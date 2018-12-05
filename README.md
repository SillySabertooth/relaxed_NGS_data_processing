# relaxed_NGS_data_processing
This bash-script (with relaxed processing of NGS data by simple programs) was made up for russian speaking people that new with NGS.

Sorry, but i have writen its for russian speaking and don't want to translate on Eng. The code is simple - all is clear.  
If somebody interest in the script (i have no idea why it can be),  
i will answer to all questions immediately - danatyermakovich@gmail.com<br/>
And maybe i might put some English explanation inside the script.

####################RU######################

Релаксированный алгоритм обработки данных NGS  
от fastq до vcf при помощи простых команд


1. Для работы вам потребуются следующие программы определенных версий  
установленные на операционной системе GNU/Linux  
(наверное, где возможно, лучше скачивать binary и т.п., ибо source нужно компилировать  
и у вас может что-то пойти не так)

Trimmomatic-0.36  
http://www.usadellab.org/cms/?page=trimmomatic

bowtie2-2.3.3 #linux-x84_64  
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3/

Samtools 1.5 & Bcftools 1.5  
https://sourceforge.net/projects/samtools/files/samtools/1.5/

vcffilter #vcflib  
https://github.com/vcflib/vcflib

ANNOVAR  
http://annovar.openbioinformatics.org/en/latest/user-guide/download/  
#annovar досутпен только последней версии  
#следует проверить команду в скрипте, вдруг они добавили изменения  
#у меня была устновлена версия 2016-02-01

на компьютере должна быть установлена java

В случае использование этих программ других версий  
в некоторых случаях опции и даже команды могут не совпадать  
(например, bcftools 0.1.19 и 1.5 вызвавают варианты разными командами - view и call соответсвенно)  
если у вас установлены другие версии  
проверьте правильность вызова команд, измените пути к программам


2. Программы следует скачать и сложить в одну папку 

я не буду заморачиваться и просить прописать путь для каждой из программ <br/>
на случай того, что они находятся в разных местах<br/>
это усложнит этот маленький скрипт и потребует много взаимодействия между<br/>
пользователем и программой<br/>
если вам тяжело сложить все программы в одну папку... 

Если вы понимаете, что здесь нет ничего сложного,<br/>
и у вас программы разбросаны по папкам и часть из них добавлена в PATH:<br/>
измените скрипт - удалите пути для "PATH"-программ<br/>
допишите конкретные пути к другим<br/>
(.jar не может быть вызван при добавлении пути его расположения в PATH)<br/>
удалите уже ненужную перменную<br/>
будет проще пользоваться, да<br/>


3. #for_newbie<br/>
Для запуска скрипта, стоит зайти в содержащую его папку через терминал<br/>
при помощи команды cd:<br/> 
например, cd /media/illumina<br/> 
#пусть он у вас лежит в этойпапке#<br/>

и затем вызвать его напрямую bash ./NGS_pipeline.sh<br/>
тут же его можно сделать исполняемым при помощи следующей команды:<br/>
chmod +x ./NGS_pipeline.sh<br/>
и запускать без слова bash -> ./NGS_pipeline.sh<br/>

также можно вызвать из любого места, указав путь сразу<br/>
bash /media/illumina/NGS_pipeline.sh<br/>
если скрипт уже исполняемый -> /media/illumina/NGS_pipeline.sh<br/>

При необходимости изменения опций программ и иного редактирования<br/>
можно открыть файл в любом текстовом редакторе<br/>
(желательно с подсветкой bash синтаксиса; тот же gedit на GNU/Linux).<br/>


4. Я не обрабатывал исключения - это слишком долго и бессмысленно<br/>
весь stderr относительно запуска скрипта <br/>
будет постпупать в терминал<br/>
и вы можете сами догадаться, почитав<br/>
где вы забыли пробел, где путь не существует, где папка существует и т.д.<br/>

почти все используемые программы записывают файлы в stdout<br/>
а лог - в stderr; чтобы на засорять терминал и для его (лога) сохранности<br/>
stderr перенаправлен в лог-файлы<br/>
и вы можете всегда почитать, как был выполнен тот или иной этап<br/>
правильно/не_правильно и с статисткиой рассчетов <br/>


5. В начале алогритма есть маленькая команда<br/>
превращающая согласно регулярному выражению стандартные названия файлов<br/>
поступающие с прибора MiSeq Illumina в имена необходимой формы<br/>
если имена имеют другую стуктуру - в надлежащий вид их можно превратить<br/>
просто переименовав в ручную или изменив регулярное выражение внутри скрипта<br/>


6. Коротко по программам:

А) Trimmomatic обрезает адаптеры на основании информации из файла TruSeq3-PE.fa<br/>
Остальные опции программы усредненнны.<br/>
Если использовалась другая пробоподготовка и другие адаптеры<br/>
можно изменить файл - найти более подходящий в папке Trimmomatic-0.36/adapters <br/>
но мне кажется, это не слишком существенно,<br/>
чтобы прописывать отдельно ввод файла с адаптерами, учитывая что не всегда очевидно, какой использовать<br/>
и огромной разницы при варьировании файлов получено не было<br/>

В) Новый референс - тогда с ним идет следующуая работа: samtools делает почти мгновенную индексацию, в то время как bowtie2-build долго строит индексы для выравнивания. Индексация референса для выравнивания занимает много времени и в этот момент онлайн-отгрузка информации в файлы может подвисать - может не отмечаться увеличения файлов-индексов в Krusader. Но, скорее всего, программа не зависла. У меня индексация рефернса hg19 размером в 3 Gb заняла 2 часа на одном ядре в 2 Gb. Подвисала раза 3. B процессе можно анализировать Log.txt.
 
C) Bowti2 выравнивает очищенные прочтения. Алгоритмы выравнивания парных и непарных прочтений разнятся между собой. Так как после trimmomatic получается и те, и другие, их следует выравнивать двумя разными вызовами bowtie2.<br/> 
Весь вопрос состоит в том, насколько это целесообразно делать. Количество непарных прочтений обычно невелико. Парные и непарные прочтения выравниваются различными алгоритмами. Соответсвенно, при сливании получаемых из этих алгоритмов файлов в один может быть привнесена некоторая погрешность. Объединение файлов требуется для дальнейшего анализа.<br/>
Опять таки, вопрос в том, насколько вы хотите выжать максимум из ваших данных и насколько это соотносится с вносимой ошибкой. В случае чего, вы можете поробовать оба варианта

D) samtools конвертирует и вызывает статистику. Здесь происходит этап удаления дуплексов. Он необходим, но из-за необоснованных опасений потерять данные я не производил его раньше на своих данных. Предлагаю обработать данные с/без удаления, сравнить результаты по кол-ву и качеству выхода и решить для себя. Думаю в дальнейшем я буду удалять дуплексы. Но вам следует определиться самостоятельно.

E) Samtools собирает позиции, Bcftools вызывает варианты. В данном месте происходит фильтрация по глубине покрытия > 15, качества картирования > 30 и уверененности прибора > 40:<br/>
вот команда (строка 245)<br/>
vcffilter -f "DP > 15 & MQ > 30 & QUAL > 40 " ${i}.raw.vcf > $Path/${i}.filt.vcf<br/>
если хотите изменить параметр - откройте скрипт и измените цифры<br/>
также можно вызвать справку по vcffilter <br/>
(путем указания команды)<br/>
и узнать, по каким параметрам еще можно фильтровать и как удалять лишние<br/>

F) annovar аннотирует :)<br/>
Перед запуском скрипта НЕОБХОДИМО глянуть мануал ANNOVARa, открыть скрипт, сравнить команды. В моей версии ANNOVAR аннотирует при помощи скаченных баз данных - их нужно предварительно все скачать!! Нужно внимательно проверить опции в мануале и скрипте, чтобы не было никаких ляпов.

В команде после -protocol идут названия баз данных<br/>
после чего в -operation идут опции<br/>
кол-во баз данных должно равняться опциям!<br/>

скаченные базы данных следует распаковать в папку annovar/humandb/<br/>
команда для скачивания бд выглядит примерно так <br/>
(стоит уточнить на сайте; .pl - скрипт, лежащий в папке annovara)<br/>
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar <база_данных> humandb/<br/>
выбрать базы данных можно на сайте<br/>
https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/download/#-for-filter-based-annotation

если вы не знаете, какие скачать - скачивайте и распаковывайте те, которые находятся у нас в скрипте<br/>
refGene,esp6500siv2_all,avsnp147,clinvar_20170501,revel,intervar_20170202,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,gnomad_genome

