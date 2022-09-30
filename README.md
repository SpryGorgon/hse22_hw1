# ВШЭ майнор "Биоинформатика" 2022/23 2 год. ДЗ 1

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

### 3. Получение статистики по исходным чтениям при помощи FastQC и MultiQC
mkdir fastqc

ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}


mkdir multiqc

multiqc -o multiqc fastqc

*На картинках pe означает paired-end, mp означает mate-pairs*
![image](https://user-images.githubusercontent.com/52814490/193300648-925aec9c-5e91-4fc6-87ec-f3a024b680ef.png)
![image](https://user-images.githubusercontent.com/52814490/193301179-5910d7e8-03dc-46e3-85b7-61cd71628143.png)
![image](https://user-images.githubusercontent.com/52814490/193305653-a38bd4fe-8d75-411b-8bd7-243cfff9d10d.png)
![image](https://user-images.githubusercontent.com/52814490/193305702-bc9c0802-7034-4ae3-963c-5cdac5bad701.png)
![image](https://user-images.githubusercontent.com/52814490/193305751-7305841a-4b23-48fa-9e60-544227d27dc6.png)
![image](https://user-images.githubusercontent.com/52814490/193305849-065ed74a-77fa-4c5d-a352-1d3e11a7d859.png)
![image](https://user-images.githubusercontent.com/52814490/193305894-6c77913f-47cb-49e8-8746-44cc53a77e6e.png)
![image](https://user-images.githubusercontent.com/52814490/193306429-58eb5572-bf69-43bd-b523-b14703515d75.png)
![image](https://user-images.githubusercontent.com/52814490/193306013-f59018d2-add0-4060-ba4a-23ad13e392dd.png)
![image](https://user-images.githubusercontent.com/52814490/193306490-64d192fe-2934-4821-878d-562efef95349.png)
![image](https://user-images.githubusercontent.com/52814490/193306207-3ad662b4-1cec-4d17-93e3-9562357658a2.png)
![image](https://user-images.githubusercontent.com/52814490/193306266-3486a67f-3ebe-4b24-94f3-47812cb1d049.png)

### 4. Подрезание чтений по качеству, удаление адаптеров
platanus_trim sub_pe_1.fastq sub_pe_2.fastq

platanus_internal_trim sub_mp_1.fastq sub_mp_2.fastq


rm sub_pe_1.fastq sub_pe_2.fastq sub_mp_1.fastq sub_mp_2.fastq

mv sub_pe_1.fastq.trimmed sub_pe_1.fastq

mv sub_pe_2.fastq.trimmed sub_pe_2.fastq

mv sub_mp_1.fastq.int_trimmed sub_mp_1.fastq

mv sub_mp_2.fastq.int_trimmed sub_mp_2.fastq

### 5. Получение статистики подрезанных чтений при помощи FastQC и MultiQC
mkdir fastqc_trimmed

ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc_trimmed {}


mkdir multiqc_int_trimmed

multiqc -o multiqc_int_trimmed fastqc_trimmed

*На картинках pe означает paired-end, mp означает mate-pairs*
![image](https://user-images.githubusercontent.com/52814490/193307856-cb98ab3b-e2d5-4233-b709-58ed76eab0d9.png)
![image](https://user-images.githubusercontent.com/52814490/193307907-225b7d57-bec0-4d9d-bb2b-9a469b5903a9.png)
![image](https://user-images.githubusercontent.com/52814490/193307946-bbd5ca7e-f0bf-4143-86c9-2856bffdcaf4.png)
![image](https://user-images.githubusercontent.com/52814490/193307979-acaa9d63-95bb-4581-a503-32488f9ff50f.png)
![image](https://user-images.githubusercontent.com/52814490/193308017-c02cd712-e81b-437e-aaca-e8a23d0090bb.png)
![image](https://user-images.githubusercontent.com/52814490/193308086-30b0672d-f312-4491-9eb3-b4c2f0eb1c66.png)
![image](https://user-images.githubusercontent.com/52814490/193308134-4144245c-1399-435b-b471-1c1ba12eeb72.png)
![image](https://user-images.githubusercontent.com/52814490/193308203-935a1a81-634f-4feb-a41a-1187bdb13843.png)
![image](https://user-images.githubusercontent.com/52814490/193308247-6688a2cb-8665-49b3-b7b1-4d73d8e13493.png)
![image](https://user-images.githubusercontent.com/52814490/193308285-680eac0f-fda6-4f08-9fbf-fa41486054b3.png)
![image](https://user-images.githubusercontent.com/52814490/193308340-19eba07e-8496-4f09-b32b-29f43cd43ee8.png)
![image](https://user-images.githubusercontent.com/52814490/193308382-372aca23-a9e5-4fa8-8604-7eb857a22fbb.png)

### 6. Сборка контигов из подрезанных чтений
platanus assemble -o pe -f sub_pe_[12].fastq -t 16 -m 56 2> assemble_pe.log

Jupyter notebook с анализом полученных контигов: src\Основное задание\contig_analysis.ipynb

### 7. Сборка скаффолдов из контигов и подрезанных чтений
platanus scaffold -o scaffold -c pe_contig.fa -b pe_contigBubble.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> scaffold.log

Jupyter notebook с анализом полученных скаффолдов: src\Основное задание\scaffold_analysis.ipynb

### 8. Уменьшение количества гэпов с помощью подрезанных чтений
platanus gap_close -o gap_close -c scaffold_scaffold.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> gap_close.log


## Бонус
Для бонусного задания возьмём 500 тысяч чтений для paired-end и 150 тысяч чтений для mate-pairs

### 1. Создание символьных ссылок на файлы с исходными данными
ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq

ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq

ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq

ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq

### 2. Выбираем случайно 500 тысяч чтений типа paired-end и 150 тысяч чтений типа mate-pairs (seed=1218)
seqtk sample -s 1218 oil_R1.fastq 500000 > sub_pe_1.fastq

seqtk sample -s 1218 oil_R2.fastq 500000 > sub_pe_2.fastq

seqtk sample -s 1218 oilMP_S4_L001_R1_001.fastq 150000 > sub_mp_1.fastq

seqtk sample -s 1218 oilMP_S4_L001_R2_001.fastq 150000 > sub_mp_2.fastq

### 3. Получение статистики по исходным чтениям при помощи FastQC и MultiQC
mkdir fastqc

ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}


mkdir multiqc

multiqc -o multiqc fastqc

*На картинках pe означает paired-end, mp означает mate-pairs*
![image](https://user-images.githubusercontent.com/52814490/193309565-f54d9819-9b05-4bca-9683-21a2eadfb233.png)
![image](https://user-images.githubusercontent.com/52814490/193309629-350d181e-e8eb-4658-8bad-7ec9f40d9e28.png)
![image](https://user-images.githubusercontent.com/52814490/193309668-00f2ad29-4919-4730-8489-c3ff290a390e.png)
![image](https://user-images.githubusercontent.com/52814490/193309716-f852af8f-56ee-427d-8fe7-c8725ed9e50e.png)
![image](https://user-images.githubusercontent.com/52814490/193309761-ba618e41-6a36-4813-aedd-ff59d03f8515.png)
![image](https://user-images.githubusercontent.com/52814490/193309803-3265497a-f493-4cd6-86a1-72d3a3dcdffe.png)
![image](https://user-images.githubusercontent.com/52814490/193309837-a097e534-12d0-4fe8-bb93-c2433495061b.png)
![image](https://user-images.githubusercontent.com/52814490/193309856-2ed734f2-47e3-4c97-b678-d9e2a2e41fa3.png)
![image](https://user-images.githubusercontent.com/52814490/193309897-4ceffbe6-d66f-4237-9831-c46d1385ddbb.png)
![image](https://user-images.githubusercontent.com/52814490/193309922-20ab9c07-1656-47ad-bf99-76d469aca1f5.png)
![image](https://user-images.githubusercontent.com/52814490/193309967-ad9435c6-6afb-4ae0-913a-a652108af265.png)
![image](https://user-images.githubusercontent.com/52814490/193310000-e93aae59-ac95-4574-9ae5-a646c73ec190.png)

### 4. Подрезание чтений по качеству, удаление адаптеров
platanus_trim sub_pe_1.fastq sub_pe_2.fastq

platanus_internal_trim sub_mp_1.fastq sub_mp_2.fastq


rm sub_pe_1.fastq sub_pe_2.fastq sub_mp_1.fastq sub_mp_2.fastq

mv sub_pe_1.fastq.trimmed sub_pe_1.fastq

mv sub_pe_2.fastq.trimmed sub_pe_2.fastq

mv sub_mp_1.fastq.int_trimmed sub_mp_1.fastq

mv sub_mp_2.fastq.int_trimmed sub_mp_2.fastq

### 5. Получение статистики подрезанных чтений при помощи FastQC и MultiQC
mkdir fastqc_trimmed

ls sub*.fastq | xargs -P 4 -tI{} fastqc -o fastqc_trimmed {}


mkdir multiqc_int_trimmed

multiqc -o multiqc_int_trimmed fastqc_trimmed

*На картинках pe означает paired-end, mp означает mate-pairs*
![image](https://user-images.githubusercontent.com/52814490/193310119-cb9f64f7-567e-4f41-90f2-b407f44b1d4d.png)
![image](https://user-images.githubusercontent.com/52814490/193310187-5e7f8c33-efe9-40a3-8282-6afc09f544ea.png)
![image](https://user-images.githubusercontent.com/52814490/193310225-a043238a-810f-48a1-9efd-e49e91171eab.png)
![image](https://user-images.githubusercontent.com/52814490/193310278-9cce9361-c1fc-4334-afc0-a37d89a8ba3d.png)
![image](https://user-images.githubusercontent.com/52814490/193310312-82333113-3f75-4c38-8aae-44e36e2a120f.png)
![image](https://user-images.githubusercontent.com/52814490/193310357-1bf035ff-b3dc-4826-b227-57325d70987a.png)
![image](https://user-images.githubusercontent.com/52814490/193310403-2171049e-2e2b-4803-a725-8bf2fb7ee67a.png)
![image](https://user-images.githubusercontent.com/52814490/193310441-1930917b-c435-4214-a93a-0c2c79438362.png)
![image](https://user-images.githubusercontent.com/52814490/193310490-ffe90561-a623-4d9d-a0ba-2d9a5ad1284d.png)
![image](https://user-images.githubusercontent.com/52814490/193310515-3fee368a-d322-4d5c-a907-111d3f49857e.png)
![image](https://user-images.githubusercontent.com/52814490/193310556-cda0a776-d506-4a57-81b2-8139121b3e41.png)
![image](https://user-images.githubusercontent.com/52814490/193310581-23ade645-c091-4bae-ac73-bdb058b1cf25.png)

### 6. Сборка контигов из подрезанных чтений
platanus assemble -o pe -f sub_pe_[12].fastq -t 16 -m 56 2> assemble_pe.log

Jupyter notebook с анализом полученных контигов: src\Бонусное задание\contig_analysis.ipynb

### 7. Сборка скаффолдов из контигов и подрезанных чтений
platanus scaffold -o scaffold -c pe_contig.fa -b pe_contigBubble.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> scaffold.log

Jupyter notebook с анализом полученных скаффолдов: src\Бонусное задание\scaffold_analysis.ipynb

### 8. Уменьшение количества гэпов с помощью подрезанных чтений
platanus gap_close -o gap_close -c scaffold_scaffold.fa -IP1 sub_pe_1.fastq sub_pe_2.fastq -OP2 sub_mp_1.fastq sub_mp_2.fastq -t 16 2> gap_close.log

### 9. Сравнение сборки, полученной из меньшего количества чтений, со сборкой из основной части
Общее количество контигов: было 602, стало 3462

Общее количество скаффолдов: было 68, стало 490

Значение N50 для контигов: было 299, стало 1735

Значение N50 для скаффолдов: было 1, стало 4

Длина самого длинного контига: было 143'791, стало 27'177

Длина самого длинного скаффолда: было 3'885'015, стало 887'110

Количество гэпов в самом длинном скаффолде: было 146, стало 641

Общая длина гэпов в самом длинном скаффолде: было 6581, стало 20352

Количество гэпов в самом длинном скаффолде после gap-close: было 23, стало 146

Общая длина гэпов в самом длинном скаффолде после gap-close: было 1246, стало 9052
