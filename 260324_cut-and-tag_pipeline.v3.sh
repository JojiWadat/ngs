#260323_cut-and-tag_pipeline.sh
#260324_cut-and-tag_pipeline.v2.sh-領域別のピーク分布図における縦軸を「ピークの割合(%)」に変更した。
#概要は末尾に記載
#!/bin/bash
set -euo pipefail
# server: horn

# =========================================================
# === Parameters manually configured ===
# =========================================================
test=0 # yes: 1 or no: 0
project="260220_VH01729_163_MH_DB" #プロジェクト名を指定（例: 260220_VH01729_163_MH_DB）。この名前をもとにディレクトリ構造が作られます。
seq="PE" #シーケンスの種類を指定。ペアエンドの場合は "PE"、シングルエンドの場合は "SE" を入力してください。これにより、FastQCのリード数の計算方法が変わります。 
threads=30 #使用するスレッド数を指定。マッピングやその他の計算で並列処理を行う際に使用されます。サーバーのCPUコア数に合わせて適切な値を設定してください。   
build="mm39_dm6" #hg38 or mm10 or mm39 #リファレンスゲノムビルドを指定。Bowtie2のインデックスとブラックリストもこのビルドに合わせて選択されます。
trimGalore_params="--nextseq 20 --cores 7"

# === Directories ===
base="/home/wada/project/${project}" #プロジェクトのベースディレクトリ。ここに "rawdata", "trim", "bam", "bw_spikein", "spikein-ratio", "bedgraph", "seacr", "figure" などのサブディレクトリが作成されます。
export BASE_DIR="$base"
raw_dir="${base}/rawdata"
txt_file="${base}/sample_info/sample_info.txt" #サンプル情報を記載したテキストファイルのパス。タブ区切りで、1列目にR1のファイル名、2列目にR2のファイル名、3列目にサンプルのプレフィックス（例: MH_P2_K4_rep1）を記載してください。ファイル名は "rawdata" ディレクトリ内のものと一致させてください。 このファイルをもとに、FastQCの出力やマッピング結果のファイル名が決まります。

# === src ===
sing="singularity exec --bind /work,/home"
bowtie2_sif="/work/SingularityImages/churros.1.2.2.sif"
BL="/work/rtao/blacklist/mm39.excluderanges.bed" #ブラックリストファイルのパス。Bowtie2でマッピングした後、samtools view -b -L "$BL" -o output.rmBL.bam input.bam でブラックリスト領域を除去します。ビルドに合わせたブラックリストを使用してください（例: mm39の場合は mm39.excluderanges.bed）。   
bw2index="/home/wada/database/bowtie2-indexes/UCSC-${build}" #Bowtie2のインデックスファイルのパス。リファレンスゲノムビルドに合わせたインデックスを指定してください（例: hg38の場合は UCSC-hg38）。このディレクトリには、bowtie2-buildで作成されたインデックスファイル（.1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2）が含まれている必要があります。  
bowtie2_params="-p $threads -N 1 -X 2000 --no-mixed --no-discordant -x $bw2index"
bowtie2_suffix="bowtie2-${build}-N1-X2000-noMixed-noDisc"
bamCoverage_param="-bs 1 -e -p $threads" # Spike-inを用いる場合は以下のようにすると2重正規化になるので注意してください。
#bamCoverage_param="--normalizeUsing BPM -bs 1 -e -p $threads"
suffix_bw=bamCov-bs1-extend-BPM.bigwig

# Setup directories
qc_1st_dir="$base/QC/FastQC/raw"; mkdir -p "$qc_1st_dir"
qc_2nd_dir="$base/QC/FastQC/trim"; mkdir -p "$qc_2nd_dir"
trim_dir="$base/trim"; mkdir -p "$trim_dir"
bam_dir="$base/bam"; mkdir -p "$bam_dir"
bw_spikein="$base/bw_spikein"; mkdir -p "$bw_spikein"
ratio_dir="$base/spikein-ratio"; mkdir -p "$ratio_dir"
bedgraph_dir="$base/bedgraph"; mkdir -p "$bedgraph_dir"
seacr_dir="$base/seacr"; mkdir -p "$seacr_dir"
fig_dir="$base/figure"; mkdir -p "$fig_dir"

# =========================================================
# === Functions ===
# =========================================================
get_bname() { bname=$(basename "$1"); }
get_dirname() { dirname=$(dirname "$1"); }
echo_message() { echo -e "$(LANG=C date '+%Y-%m-%d %H:%M:%S'):\t$*"; }
check_and_execute() {
    local file="$1"
    local command_local="$2"
    if [[ "$test" == 1 ]]; then
        echo_message "$(echo "$command_local" | awk '{print $1}'):\tThis is test running"
    else
        if [[ ! -e "$file" ]]; then
            eval "$command_local"
            echo -e ""
        else
            echo_message "file '$file' already exists. Skip this step."
            echo -e ""
        fi
    fi
}
run_fastqc() { fastqc --nogroup -t "$3" -o "$2" "$1" 2>/dev/null; }
run_trim_galore() { trim_galore ${trimGalore_params} -o "$3" --basename "$4" --paired "$1" "$2"; }
run_bowtie2() { ${sing} ${bowtie2_sif} bowtie2 ${bowtie2_params} -1 "$1" -2 "$2" 2> "$4" | samtools view -bS -@ "$5" -h - | samtools sort -@ "$5" -o "$3" -; }
remove_chrM() { samtools view -@ "$threads" -h "$1" | awk 'BEGIN{OFS="\t"} /^@/ {print; next} $3!="chrM"{print}' | samtools sort -@ "$threads" -o "$2" -; }
index_bam_if_needed() { if [[ ! -e "${1}.bai" && ! -e "${1%.bam}.bai" ]]; then samtools index -@ "$threads" "$1"; fi; }
Filtering_MAPQ30_indexing() { samtools view -@ "$threads" -b -q 30 $1 > $2; index_bam_if_needed "$2"; }
remove_PCRdup() {
    local odir="$(dirname "$2")"
    samtools collate -@ "$threads" -o "$odir/namecollate.bam" "$1"
    samtools fixmate -@ "$threads" -m "$odir/namecollate.bam" "$odir/fixmate.bam"
    samtools sort -@ "$threads" -o "$odir/positionsort.bam" "$odir/fixmate.bam"
    samtools markdup -r -S -@ "$threads" "$odir/positionsort.bam" "$2"
    rm "$odir/namecollate.bam" "$odir/fixmate.bam" "$odir/positionsort.bam"
}
remove_blacklist() { bedtools intersect -v -abam "$1" -b "$BL" -ubam > "$2"; index_bam_if_needed "$2"; }
obtain_read_fromFastQC() {
    local fastqc_data="${1%.zip}/fastqc_data.txt"
    unzip -q "$1" -d "$(dirname "$1")"
    n=$(awk -F "\t" 'NR==7 {print $2}' "$fastqc_data")
    if [[ "$seq" == "PE" ]]; then n=$((n * 2)); fi
    n_cs=$(printf "%'d\n" "$n")
    rm -r "${1%.zip}"
}
count_percent() { echo "$1 $2" | awk '{printf "%.1f\n", 100*$1/$2}'; }
convert_to_bw() { bamCoverage --bam "$1" -o "$2" $bamCoverage_param; }

# =========================================================
# === 1. Processing Pipeline (FastQC to BigWig) ===
# =========================================================
declare -a FNAME_R1_RAW FNAME_R2_RAW PREFIX MAPPED MAPQ30 rmCHRM rmPCRDUP rmBL
while read -r line; do
    read -ra ARY <<< "$line"
    FNAME_R1_RAW+=(${ARY[0]})
    FNAME_R2_RAW+=(${ARY[1]})
    PREFIX+=(${ARY[2]})
done < "$txt_file"

echo_message "START\t: FastQC and trimming"
for i in "${!PREFIX[@]}"; do
    # FastQC
    check_and_execute "$qc_1st_dir/${FNAME_R1_RAW[i]%.fastq.gz}_fastqc.html" "run_fastqc '$raw_dir/${FNAME_R1_RAW[i]}' '$qc_1st_dir' '$threads'"
    check_and_execute "$qc_1st_dir/${FNAME_R2_RAW[i]%.fastq.gz}_fastqc.html" "run_fastqc '$raw_dir/${FNAME_R2_RAW[i]}' '$qc_1st_dir' '$threads'"
    
    # Trim Galore
    trim_r1="$trim_dir/${PREFIX[i]}_val_1.fq.gz"
    trim_r2="$trim_dir/${PREFIX[i]}_val_2.fq.gz"
    check_and_execute "$trim_r1" "run_trim_galore '$raw_dir/${FNAME_R1_RAW[i]}' '$raw_dir/${FNAME_R2_RAW[i]}' '$trim_dir' '${PREFIX[i]}'"
    
    # FastQC trimmed
    check_and_execute "$qc_2nd_dir/${PREFIX[i]}_val_1_fastqc.html" "run_fastqc '$trim_r1' '$qc_2nd_dir' '$threads'"
    check_and_execute "$qc_2nd_dir/${PREFIX[i]}_val_2_fastqc.html" "run_fastqc '$trim_r2' '$qc_2nd_dir' '$threads'"

    # Mapping & Filtering
    mapping_out="$bam_dir/${PREFIX[i]}.sort.bam"
    check_and_execute "$mapping_out" "run_bowtie2 '$trim_r1' '$trim_r2' '$mapping_out' '$bam_dir/${PREFIX[i]}-${bowtie2_suffix}.bowtie2Log.txt' '$threads'"
    
    mapq30_out="${mapping_out%.sort.bam}.MAPQ30.sort.bam"
    check_and_execute "$mapq30_out" "Filtering_MAPQ30_indexing '$mapping_out' '$mapq30_out'"
    
    rmchrm_out="${mapq30_out%.sort.bam}.rmChrM.sort.bam"
    check_and_execute "$rmchrm_out" "remove_chrM '$mapq30_out' '$rmchrm_out'"
    
    rmpcr_out="${rmchrm_out%.sort.bam}.rmPCRdup.sort.bam"
    check_and_execute "$rmpcr_out" "remove_PCRdup '$rmchrm_out' '$rmpcr_out'"
    
    rmbl_out="${rmpcr_out%.sort.bam}.rmBL.sort.bam"
    check_and_execute "$rmbl_out" "remove_blacklist '$rmpcr_out' '$rmbl_out'"
    rmBL+=("$rmbl_out")
done

# =========================================================
# === 2. Spike-in Normalization ===
# =========================================================
echo_message "START\t: Spike-in normalization"
ratio_file="$ratio_dir/spikein_ratio.tsv"
echo -e "sample\tmouse_reads\tfly_reads\tmouse_ratio\tfly_ratio" > "${ratio_file}"
for bam in "${bam_dir}"/*.rmBL.sort.bam; do
    name=$(basename "$bam" .rmBL.sort.bam)
    dm6=$(samtools idxstats "$bam" | awk '$1 ~ /^fly_/ {sum+=$3} END{print sum+0}')
    mm39=$(samtools idxstats "$bam" | awk '$1 !~ /^fly_/ {sum+=$3} END{print sum+0}')
    total=$((dm6 + mm39))
    if [[ $total -gt 0 ]]; then
        mouse_ratio=$(awk -v m=$mm39 -v t=$total 'BEGIN{printf "%.6f", m/t}')
        fly_ratio=$(awk -v d=$dm6 -v t=$total 'BEGIN{printf "%.6f", d/t}')
        echo -e "$name\t$mm39\t$dm6\t$mouse_ratio\t$fly_ratio" >> "${ratio_file}"
    fi
done

max_dm6=$(awk 'NR>1 {if($3>max) max=$3} END{print max+0}' "$ratio_file")
if [[ "$max_dm6" -gt 0 ]]; then
    for bam in "$bam_dir"/*.rmBL.sort.bam; do
        name=$(basename "$bam" .rmBL.sort.bam)
        dm6=$(samtools idxstats "$bam" | awk '$1 ~ /^fly_/ {sum+=$3} END{print sum+0}')
        if [[ "$dm6" -gt 0 ]]; then
            scale=$(awk -v ref=$max_dm6 -v d=$dm6 'BEGIN{printf "%.10f", ref/d}')
            out_bw="$bw_spikein/${name}.spikein.bigwig"
            check_and_execute "$out_bw" "bamCoverage --bam '$bam' -o '$out_bw' $bamCoverage_param --scaleFactor '$scale'"
        fi
    done
fi

# =========================================================
# === 3. Bedgraph conversion & SEACR Peak Calling ===
# =========================================================
echo_message "START\t: BedGraph & SEACR Peak Calling"
stages=("P2" "P3" "P4" "P5" "P7")
reps=("rep1" "rep2")

for stage in "${stages[@]}"; do
    # 3-1: 対象のMarkを設定（P2のみK4とK27、他はK27のみ）
    if [[ "$stage" == "P2" ]]; then
        marks=("K4" "K27")
    else
        marks=("K27")
    fi

    for rep in "${reps[@]}"; do
        # 3-2: IgGをBedgraphに変換
        igg_bam="$bam_dir/MH_${stage}_IgG_${rep}.MAPQ30.rmChrM.rmPCRdup.rmBL.sort.bam"
        igg_bg="$bedgraph_dir/MH_${stage}_IgG_${rep}.bedgraph"
        if [[ -f "$igg_bam" && ! -f "$igg_bg" ]]; then
            echo_message "Converting $igg_bam to bedgraph..."
            bedtools genomecov -bg -ibam "$igg_bam" > "$igg_bg"
        fi

        # 3-3: ターゲット(K4/K27)をBedgraphに変換し、SEACRを実行
        for mark in "${marks[@]}"; do
            tgt_bam="$bam_dir/MH_${stage}_${mark}_${rep}.MAPQ30.rmChrM.rmPCRdup.rmBL.sort.bam"
            tgt_bg="$bedgraph_dir/MH_${stage}_${mark}_${rep}.bedgraph"
            
            if [[ -f "$tgt_bam" && ! -f "$tgt_bg" ]]; then
                echo_message "Converting $tgt_bam to bedgraph..."
                bedtools genomecov -bg -ibam "$tgt_bam" > "$tgt_bg"
            fi
            
            seacr_out_prefix="$seacr_dir/MH_${stage}_${mark}_${rep}_seacr"
            if [[ -f "$tgt_bg" && -f "$igg_bg" && ! -f "${seacr_out_prefix}.stringent.bed" ]]; then
                echo_message "Running SEACR for ${stage}_${mark}_${rep}..."
                SEACR_1.3.sh "$tgt_bg" "$igg_bg" non stringent "$seacr_out_prefix"
            fi
        done
    done
done

# =========================================================
# === 4 & 5. R Processing (Setup & Analysis) ===
# =========================================================
echo_message "START\t: R Analysis Pipeline"
R_SCRIPT_PATH="$base/run_seacr_analysis.R"

cat << 'EOF' > "$R_SCRIPT_PATH"
# --- R パッケージのインストールと読み込み ---
cran_packages <- c("ggplot2", "dplyr", "eulerr", "BiocManager")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos="http://cran.us.r-project.org")
}
bioc_packages <- c("ChIPseeker", "TxDb.Mmusculus.UCSC.mm39.refGene", "GenomicRanges")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE)
}

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(eulerr)

# --- 環境設定 ---
base_dir <- Sys.getenv("BASE_DIR")
seacr_dir <- file.path(base_dir, "seacr")
fig_dir <- file.path(base_dir, "figure")
setwd(seacr_dir)

txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
mouse_chrs <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")

# 解析対象の組み合わせリスト
targets <- c("P2_K4", "P2_K27", "P3_K27", "P4_K27", "P5_K27", "P7_K27")

for (tgt in targets) {
  rep1_file <- paste0("MH_", tgt, "_rep1_seacr.stringent.bed")
  rep2_file <- paste0("MH_", tgt, "_rep2_seacr.stringent.bed")
  
  if (file.exists(rep1_file) && file.exists(rep2_file)) {
    message("Processing: ", tgt)
    
    # ピーク読み込みと染色体フィルタリング
    peak_rep1 <- readPeakFile(rep1_file)
    peak_rep2 <- readPeakFile(rep2_file)
    peak_rep1 <- peak_rep1[seqnames(peak_rep1) %in% mouse_chrs]
    peak_rep2 <- peak_rep2[seqnames(peak_rep2) %in% mouse_chrs]
    if (length(peak_rep1) == 0 || length(peak_rep2) == 0) {
        message("Skipping ", tgt, ": No mouse peaks found in one or both replicates.")
        next
    }
    
    # =========================================================
    # --- グラフ①：面積比例のベン図 (eulerr) ---
    # =========================================================
    # GRangesを用いた正確なオーバーラップ計算
    n_rep1_total <- length(peak_rep1)
    n_rep2_total <- length(peak_rep2)
    
    # Rep1とRep2の重なるピーク数を算出 (基準をRep1に置く)
    n_overlap <- length(subsetByOverlaps(peak_rep1, peak_rep2))
    
    # それぞれのユニークなピーク数を算出
    n_rep1_only <- n_rep1_total - n_overlap
    # 1対多の重なりを補正するため、Rep2のみの数は (Rep2総数 - Rep1に重なるRep2数) で計算
    n_rep2_only <- n_rep2_total - length(subsetByOverlaps(peak_rep2, peak_rep1))
    
    # Euler図のフィッティング
    fit <- euler(c("Rep1" = n_rep1_only, "Rep2" = n_rep2_only, "Rep1&Rep2" = n_overlap))
    
    venn_path <- file.path(fig_dir, paste0("MH_", tgt, "_SEACR_Overlap_Euler.pdf"))
    pdf(venn_path, width = 5, height = 5)
    p_euler <- plot(fit, 
                    quantities = TRUE, 
                    fill = c("#66C2A5", "#FC8D62"), 
                    alpha = 0.5,
                    labels = c("Rep1", "Rep2"),
                    main = paste(tgt, "Peak Overlap (Area Proportional)"))
    print(p_euler) # pdfデバイスへ明示的に出力
    dev.off()
    
    # =========================================================
    # --- グラフ②：ピークの「分布割合 (Percentage)」の棒グラフ ---
    # =========================================================
    # アノテーションの実行
    peak_anno_rep1 <- annotatePeak(peak_rep1, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
    peak_anno_rep2 <- annotatePeak(peak_rep2, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE)
    
    process_anno <- function(anno_obj, rep_name) {
      as.data.frame(anno_obj) %>%
        mutate(
          Region = case_when(
            grepl("Promoter", annotation) ~ "TSS",
            grepl("Exon|Intron", annotation) ~ "Genic",
            grepl("Downstream", annotation) ~ "TES",
            grepl("Distal Intergenic", annotation) ~ "Intergenic",
            TRUE ~ "Other"
          ),
          Replicate = rep_name
        ) %>%
        filter(Region != "Other")
    }
    
    df_rep1 <- process_anno(peak_anno_rep1, "Rep1")
    df_rep2 <- process_anno(peak_anno_rep2, "Rep2")
    
    combined_df <- bind_rows(df_rep1, df_rep2)

    # 領域ごとのカウントと「割合(%)」を算出
    plot_data <- combined_df %>%
      group_by(Replicate, Region) %>%
      summarise(Count = n(), .groups = "drop") %>%
      group_by(Replicate) %>%
      mutate(Percentage = (Count / sum(Count)) * 100)

    plot_data$Region <- factor(plot_data$Region, levels = c("Intergenic", "TSS", "Genic", "TES"))

    # ggplot2で描画 (Y軸を Percentage に)
    p_dist <- ggplot(plot_data, aes(x = Region, y = Percentage, fill = Replicate)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
      theme_classic() +
      labs(
        title = paste(tgt, "Peak Distribution by Region (SEACR)"),
        x = "Genomic Region",
        y = "Percentage of Total Peaks (%)"
      ) +
      scale_fill_manual(values = c("Rep1" = "#66C2A5", "Rep2" = "#FC8D62")) +
      theme(
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    # 保存 (ファイル名もターゲットに応じて動的に変更)
    dist_path <- file.path(fig_dir, paste0("MH_", tgt, "_SEACR_Peak_Distribution.pdf"))
    ggsave(dist_path, plot = p_dist, width = 7, height = 5)
    
  } else {
    message("Skipped ", tgt, ": Missing SEACR stringent bed files.")
  }
}
EOF

# 生成したRスクリプトの実行
Rscript "$R_SCRIPT_PATH"
echo_message "All analyses Finished!!"

#===概要===#
#1. データクレンジングとマッピング
# FastQC / Trim Galore: アダプター配列の除去と品質チェック。
# Bowtie2: マウス(mm39)とショウジョウバエ(dm6)の混合リファレンスへマッピング。
# 多段階フィルタリング:
# MAPQ30: マッピング精度の低いリードを排除。
# rmChrM: ミトコンドリア由来のノイズを除去。
# markdup -r: PCR増幅バイアス（重複）を除去。
# remove_blacklist: 異常に高いシグナルが出る既知の不良領域を除去。

# 2. スパイクイン補正（定量的解析のキモ）
# ショウジョウバエ由来のリード数（fly_reads）をカウント。
# 全サンプル中で最も効率の良かったサンプルに合わせてScale Factorを計算。
# この係数を用いてBigWigファイルを作成することで、**細胞数や実験手技のムラを排除した「定量的な比較」**が可能になります。

# 3. SEACR によるピーク検出
# CUT&Tagに最適化された SEACR を使用。
# ターゲット（K4/K27）のBedGraphと、各サンプルのIgGを直接比較し、バックグラウンドより有意に高い領域を stringent モードで抽出します。

# 4. Rによる自動統計解析・可視化
# スクリプトの後半では、抽出されたピークを使って以下の解析を自動実行します。
# ベン図: Rep1とRep2でどれくらいピークが重なっているか（再現性の確認）。
# ゲノムアノテーション: ピークを「TSS」「Genic」「TES」「Intergenic」に分類。
# 領域別ピーク分布プロット: 検出されたピークがゲノム上のどの領域にどれくらいの割合(%)で分布しているかを棒グラフ化。