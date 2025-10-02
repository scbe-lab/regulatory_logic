macs2 callpeak -f BED -t ../nucf_trim100bp/atac_vir/atac_vir_nucfree.bed --nomodel --extsize 100 --shift 45 --buffer-size 50000 -g 840213658 -p 0.001 -n atac_vir_macs2peaks --outdir ./macs2_out/
cat macs2_out/atac_vir_macs2peaks_peaks.narrowPeak | mergeBed -d 20 | sortBed | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak_"NR}' > peaks_virginia.bed

