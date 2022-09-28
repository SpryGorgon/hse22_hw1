# ВШЭ майнор "Биоинформатика" 2022/23 2 год. ДЗ 1
## Файлы

## Основное задание
### 1. Создание символьных ссылок на файлы с исходными данными
ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq

### 2. Выбираем случайно 5 миллионов чтений типа paired-end и 1.5 миллиона чтений типа mate-pairs (seed=1218)
seqtk sample -s 1218 oil_R1.fastq 5000000 > sub_pe_1.fastq
seqtk sample -s 1218 oil_R2.fastq 5000000 > sub_pe_2.fastq
seqtk sample -s 1218 oilMP_S4_L001_R1_001.fastq 1500000 > sub_mp_1.fastq
seqtk sample -s 1218 oilMP_S4_L001_R2_001.fastq 1500000 > sub_mp_2.fastq

### 3. Получение статистики по чтениям при помощи FastQC
mkdir fastqc
ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}

#### Для первого чтения
![per_base_quality](https://user-images.githubusercontent.com/52814490/192885917-3692dcd3-d5d6-48e8-a940-4473da972db9.png)
![per_tile_quality](https://user-images.githubusercontent.com/52814490/192885926-dc5b969d-ce5e-4f1c-9955-94b60df7fccf.png)
![per_sequence_quality](https://user-images.githubusercontent.com/52814490/192885925-ae60c93a-4f4c-477c-8c45-7e6cf072a1b5.png)
![per_base_sequence_content](https://user-images.githubusercontent.com/52814490/192885920-9dfdfd07-c797-40f9-a10f-00e171894b21.png)
![per_sequence_gc_content](https://user-images.githubusercontent.com/52814490/192885922-3d5b52a4-5f9e-47ce-823b-9919b236fc94.png)
![per_base_n_content](https://user-images.githubusercontent.com/52814490/192885916-49df70be-2fa3-4e8c-b577-87c829b5a968.png)
![sequence_length_distribution](https://user-images.githubusercontent.com/52814490/192885928-1c2817ba-a40c-4ece-bae6-2dec2a3a6d9d.png)
![duplication_levels](https://user-images.githubusercontent.com/52814490/192885915-b66bf5ae-dec3-4792-a674-900037b6a33d.png)
![adapter_content](https://user-images.githubusercontent.com/52814490/192885911-144b2aa0-c041-42b8-847d-dbe869b44518.png)

#### Для второго чтения
![per_base_quality](https://user-images.githubusercontent.com/52814490/192886871-009927ea-9afa-49b5-9cf4-f6642fe26319.png)
![per_tile_quality](https://user-images.githubusercontent.com/52814490/192886887-0f880a74-8248-44eb-9088-2768a3197dbe.png)
![per_sequence_quality](https://user-images.githubusercontent.com/52814490/192886882-fbba3e16-59fe-4420-9b5b-be07abe02ed1.png)
![per_base_sequence_content](https://user-images.githubusercontent.com/52814490/192886873-6bc874ef-4c10-400b-8caa-745070d93bbf.png)
![per_sequence_gc_content](https://user-images.githubusercontent.com/52814490/192886876-3727702c-0353-425a-8f11-f62faa78bbd6.png)
![per_base_n_content](https://user-images.githubusercontent.com/52814490/192886867-b3566df4-7b9e-4233-9a10-3d8a6abfd5a7.png)
![sequence_length_distribution](https://user-images.githubusercontent.com/52814490/192886894-c39230f1-45b6-4f80-8d7b-e5b3cf4b9eff.png)
![duplication_levels](https://user-images.githubusercontent.com/52814490/192886865-905bcb35-ea9c-4143-8f61-60a8b8b3eaee.png)
![adapter_content](https://user-images.githubusercontent.com/52814490/192886862-d8daaffa-1932-4585-b615-0b7604245906.png)


### 4. Получение статистики по чтениям при помощи MultiQC
mkdir multiqc
multiqc -o multiqc fastqc

#### Для первого чтения


#### Для второго чтения


platanus_trim sub_pe_1.fastq sub_pe_2.fastq
platanus_internal_trim sub_mp_1.fastq sub_mp_2.fastq

rm sub_pe_1.fastq sub_pe_2.fastq sub_mp_1.fastq sub_mp_2.fastq
mv sub_pe_1.fastq.trimmed sub_pe_1.fastq
mv sub_pe_2.fastq.trimmed sub_pe_2.fastq
mv sub_mp_1.fastq.int_trimmed sub_mp_1.fastq
mv sub_mp_2.fastq.int_trimmed sub_mp_2.fastq

mkdir fastqc_trimmed
ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc_trimmed {}

mkdir multiqc_int_trimmed
multiqc -o multiqc_int_trimmed fastqc_trimmed

platanus assemble -o pe -f sub_pe_[12].fastq -t 16 -m 56 2> assemble_pe.log
??? не используется platanus assemble -o mp -f sub_mp_[12].fastq -t 16 -m 56 2> assemble_mp.log

platanus scaffold -o scaffold -c pe_contig.fa -b pe_contigBubble.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> scaffold.log

platanus gap_close -o gap_close -c scaffold_scaffold.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> gap_close.log

Для бонуса возьмём 500к чтений для paired-end и 150к чтений для mate-pairs

ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq

seqtk sample -s 1218 oil_R1.fastq 500000 > sub_pe_1.fastq
seqtk sample -s 1218 oil_R2.fastq 500000 > sub_pe_2.fastq
seqtk sample -s 1218 oilMP_S4_L001_R1_001.fastq 150000 > sub_mp_1.fastq
seqtk sample -s 1218 oilMP_S4_L001_R2_001.fastq 150000 > sub_mp_2.fastq

mkdir fastqc
ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}

mkdir multiqc
multiqc -o multiqc fastqc

platanus_trim sub_pe_1.fastq sub_pe_2.fastq
platanus_internal_trim sub_mp_1.fastq sub_mp_2.fastq

rm sub_pe_1.fastq sub_pe_2.fastq sub_mp_1.fastq sub_mp_2.fastq
mv sub_pe_1.fastq.trimmed sub_pe_1.fastq
mv sub_pe_2.fastq.trimmed sub_pe_2.fastq
mv sub_mp_1.fastq.int_trimmed sub_mp_1.fastq
mv sub_mp_2.fastq.int_trimmed sub_mp_2.fastq

mkdir fastqc_trimmed
ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc_trimmed {}

mkdir multiqc_int_trimmed
multiqc -o multiqc_int_trimmed fastqc_trimmed

platanus assemble -o pe -f sub_pe_[12].fastq -t 16 -m 56 2> assemble_pe.log

platanus scaffold -o scaffold -c pe_contig.fa -b pe_contigBubble.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> scaffold.log

platanus gap_close -o gap_close -c scaffold_scaffold.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> gap_close.log
