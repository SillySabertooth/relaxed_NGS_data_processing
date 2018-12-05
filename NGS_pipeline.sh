#!/bin/bash

version 1.0

echo "
Релаксированный алгоритм обработки данных NGS  
советую прочитать README, там буквально страница 

Сначала укажите путь в содержащую файлы для обработки папку 
эти файлы должны иметь унифицированное название <номер_образца>.R<1или2>.fastq.gz
например, 16.R1.fastq.gz и 16.R2.fastq.gz; 
алгоритм может преобразовать названия если они имеют такой вид 
<заданный_номер>_S<номер_зад_прибором>_L001_R1_001.fastq.gz 
19_S3_L001_R1_001.fastq.gz  

если путь к файлу представляет собой /media/illumina/16.R1.fastq.gz 
то прописываемый путь дожен быть задан как /media/illumina 
путь можно скопировать любым образом и вставить при помощи Ctrl+Shift+V в терминал
"
read Path_short

echo "
Теперь введите путь к рефернсу !вместе с именем файла! 
на который будет осуществляться выравнивание  
например, media/illumina/reference/hg19.fa 
он должен быть в формате .fa, .fasta, .fna и т.п.; hg19 сборка весит ~3 Gb 
"

read Ref

echo "
И еще путь к папке со ВСЕМИ программами 
"

read Path_to_pr

echo "
Введите самый высокий номер образца в обрабатываемом сете
"

read Numb

echo "
Предпоследнее - выравнивать только парные прочтения?
Yes 
No
(если вы не знаете, что это такое - тоже вводите Yes)
"
read Answer

echo "
Последнее!!
Удалить ли дуплексы?
Yes
No
Без удаления будет больше информации. Следует помнить, что все алгоритмы имеют погрешности и может быть удаленно больше, чем следует. Лучше всего запустить программу два раза в обоих режимах и сравнить результат
+ команда удаляет вероятные дуплексы исходя из pair-ended
если вы в предыдущим поле выбрали No, то здесь лучше тоже выбрать No
наверное; результат может будет непредсказуем
"

read Rmd

#echo "
#Этим сообщением будут отделяться данные по покрытию для разных запусков скрипта
#На случай того, что вдруг вы постоянно запускаете его на одних и тех же данных в одной и той же папке
#" >> $Path_short/depth_of_final_bams.odt
#if [[ "$Rmd" == "Yes" ]]; then
#echo -e "\nИспользуемые для высчитывания BAM-файлы были получены c удалениeм дуплексов\n" >> $Path_short/depth_of_final_bams.odt
#else
#echo -e "\nИспользуемые для высчитывания BAM-файлы были получены без удаления дуплексов\n" >> $Path_short/depth_of_final_bams.odt
#fi

#####################################################
echo -e "


Этим сообщением отделяется Log для  разных запусков данного скрипта
на данных, находящихся в одной папке


" >>$Path_short/Log.txt
exec 2>> $Path_short/Log.txt
cd $Path_short 
mkdir Results/

### changing_names
for i in $(seq 1 $Numb)
do
for R in 1 2
do
rename "s/${i}_S[0-9]*_L001_R${R}_001.fastq/${i}.R${R}.fastq/" *
done
done

###bowtie_making_index_for_ref
Ref_dir=`dirname $Ref`
filename="${Ref##*/}"
filename="${filename%.*}"
if [[ -f $Ref_dir/${filename}.1.bt2 ]]; then
echo "
Ваш референс был индексирован bowtie2 раннее. Отлично!
"
else
echo "
Вы получаете это сообщение если ваш референс не индексирован
"
echo "
Начало индексации рерфернса
"
echo -e "\nBowtie2-build log to index the reference\n" >> $Path_short/Log.txt
$Path_to_pr/bowtie2-2.3.3/bowtie2-build $Ref $Ref_dir/${filename} >> $Path_short/Log.txt
echo "
Индексация закончена
"
fi

if [[ -f $Ref_dir/${filename}.fna.fai ]]; then
echo "
Также референс уже был индексирован samtools
"
else
samtools faidx $Ref
fi
#######################################

for i in $(seq 1 $Numb)
do
if [[ -f $Path_short/${i}.R1.fastq.gz ]]; then
mkdir $Path_short/${i}
Path=$Path_short/${i}
cd $Path

### trimmomatic
echo "Запускается Trimmomatic"
echo "
Отследить работу кода на каждом этапе можно по увеличению
размера файлов в рабочей папке через файловые менеджеры типа Krusader
"
echo -e "\n The Iteration for ${i}! \n\n Trim_log \n" >> $Path_short/Log.txt
java -jar $Path_to_pr/Trimmomatic-0.36/trimmomatic-0.36.jar \
PE \
$Path_short/${i}.R1.fastq.gz \
$Path_short/${i}.R2.fastq.gz \
${i}.R1_paired.fastq \
${i}.R1_unpaired.fastq \
${i}.R2_paired.fastq \
${i}.R2_unpaired.fastq \
ILLUMINACLIP:$Path_to_pr/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:45
echo "Trimmomatic окончил работу над образцом ${i}"

###Bowtie2

#bowtie2_make_paired

echo "Запускается Bowtie2"
echo -e "\nBowtie_log\n" >> $Path_short/Log.txt 

$Path_to_pr/bowtie2-2.3.3/bowtie2 \
-x $Ref_dir/${filename} \
-1 ${i}.R1_paired.fastq \
-2 ${i}.R2_paired.fastq \
-S ${i}.p.sam

if [[ "$Answer" == "No" ]]; then
echo "Запускается Bowtie2 в режиме выравнивая парных и непарных прочтений 
(инф - README п.6.C)"

#bowtie2_make_unpaired

$Path_to_pr/bowtie2-2.3.3/bowtie2 \
-x $Ref_dir/${filename} \
-U ${i}.R1_unpaired.fastq \
-U ${i}.R2_unpaired.fastq \
-S ${i}.u.sam

echo "Bowtie2 окончил обработку ${i}"

#samtools_make_BAMs
echo "Запускается Samtools"
echo -e "\nSamtools_log\n Вроде, обычно пустой в случае конвертирования файлов \n" >> $Path_short/Log.txt 

for n in p u 
do
$Path_to_pr/samtools-1.5/samtools view ${i}.${n}.sam -o ${i}.${n}.bam
$Path_to_pr/samtools-1.5/samtools sort ${i}.${n}.bam -o ${i}.${n}.s.bam
$Path_to_pr/samtools-1.5/samtools index ${i}.${n}.s.bam
done

#samtools_merge_files
$Path_to_pr/samtools-1.5/samtools merge ${i}.m.bam ${i}.p.s.bam ${i}.u.s.bam
$Path_to_pr/samtools-1.5/samtools sort ${i}.m.bam -o ${i}.m.sort.bam
$Path_to_pr/samtools-1.5/samtools index ${i}.m.sort.bam
echo "samtools создал финальный BAM-файл для ${i}"

#если вы выбрали Yes, то обработка пойдет по этому, более простому пути
else

#samtools_make_BAMs
echo "Был запущен Bowtie2 в режиме выравнивая только парных прочтений 
(инф - README п.6.C) и уже окончил обработку ${i}
Теперь запускается Samtools"
echo -e "\nSamtools_log\n Вроде, обычно пустой в случае конвертирования файлов \n" >> $Path_short/Log.txt 
$Path_to_pr/samtools-1.5/samtools view ${i}.p.sam -o ${i}.m.bam
$Path_to_pr/samtools-1.5/samtools sort ${i}.m.bam -o ${i}.m.sort.bam
$Path_to_pr/samtools-1.5/samtools index ${i}.m.sort.bam
echo "samtools создал финальный BAM-файл для ${i}"
fi

#deleting_the_duplexes
#если вы выбрали Yes, то дуплексы буду удалены
if [[ "$Rmd" == "Yes" ]]; then
echo "Происходит удаление дуплексов"
echo -e "\nSamtools_log_for_rmd\n" >> $Path_short/Log.txt 
$Path_to_pr/samtools-1.5/samtools rmdup -S ${i}.m.sort.bam ${i}.d.bam
$Path_to_pr/samtools-1.5/samtools sort ${i}.d.bam -o ${i}.m.sort.bam
$Path_to_pr/samtools-1.5/samtools index ${i}.m.sort.bam
echo "Дуплексы для ${i} удаленны"
else
echo -e "\nSamtools_log_for_rmd\n Дуплексы не удалялись \n" >> $Path_short/Log.txt 
echo "Дуплексы для ${i} не удалялись"
fi

#statistics
$Path_to_pr/samtools-1.5/samtools stats ${i}.m.sort.bam > $Path_short/Results/${i}.m_sort_bam.bc
cd $Path_short/Results/
$Path_to_pr/samtools-1.5/misc/plot-bamstats ${i}.m_sort_bam.bc -p ${i}_bamstats_plots/
mv ${i}.m_sort_bam.bc ${i}_bamstats_plots/
cd $Path/

#покрытие высчитвает отношение покрытия на весь геном
#поэтому является малоинформативным
#$Path_to_pr/samtools-1.5/samtools depth ${i}.m.sort.bam | awk '{sum+=$3} END { print "Среднее покртыие '${i}' = ",sum/NR}' >> $Path_short/depth_of_final_bams.odt

#############################################

###samtools_make_vcf
echo "Происходит вызов вариантов при помощи Samtools
сырые файлы также здесь фильтруются по DP > 15, MQ >30, QUAL > 40 
для изменения параметров фильтрации - README
"
echo -e "\nSamtools_log_for_calling\n" >> $Path_short/Log.txt 
$Path_to_pr/samtools-1.5/samtools mpileup -uvf $Ref ${i}.m.sort.bam | $Path_to_pr/bcftools-1.5/bcftools call -cv - > ${i}.raw.vcf

#vcffilter_filter_variants

$Path_to_pr/vcflib/bin/vcffilter -f "DP > 15 & MQ > 30 & QUAL > 40 " ${i}.raw.vcf > ${i}.filt_15.vcf
echo "Вызов вариантов и фильтрации закончены для образца ${i}"

###annovar_annotating :)
echo "Запускается annovar"
echo -e "\nannovar_log\n" >> $Path_short/Log.txt
echo "
Надеюсь, вы читали readme и мануал annovara, после чего 
проверили команду annovar 
для вашей версии, ваших баз данных и ваших папок
"
$Path_to_pr/annovar/table_annovar.pl ${i}.filt_15.vcf humandb/ -buildver hg19 -out ./${i}_anno -remove -protocol refGene,esp6500siv2_all,avsnp147,clinvar_20170501,revel,intervar_20170202,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,gnomad_genome -operation g,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

echo "аннотация для образца ${i} закончена"

mv ${i}.m.sort.bam  ${i}.m.sort.bam.bai ${i}.raw.vcf ${i}.filt_15.vcf ${i}_anno* $Path_short/Results/ 
cd ../
rmdir ${i}/


fi
done
