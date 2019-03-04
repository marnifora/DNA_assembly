# DNA_assembly

*English version is below*

### Narzędzie do składanie odczytów pochodzących z sekwencjonowanie DNA typu "single-end".

Polecenie do zadania znajduje się w pliku assignment2en.pdf.

#### Wersje:
Python: 3.6.6

#### Jak używać?
./assembly wejściowe_odczyty.fasta utworzone_contigi.fasta

#### Przebieg assemblacji odczytów:

1. Wczytanie readów z pliku fasta.
2. Poprawa błędnych odczytów.
    - Na podstawie zliczenia ile i jakich k-merów powstaje z danych wejściowych, wyliczane jest średnie pokrycie k-mera
        oraz odchylenie standardowe. Wartość średnia-odchylenie uznaje się za próg. Dla k-merów o pokryciu poniżej progu
         szukany jest k-mer „siostrzany” - o odległości Hamminga równej 1. Jeśli pokrycie k-meru siostrzanego jest 
         powyżej progu to wejściowy k-mer jest zamieniany w k-mer siostrzany. K-mery o pokryciu poniżej progu, dla 
         których nie istnieją odpowiednio dobre k-mery siostrzane, są zapisywane na listę wrong_kmers.
3. Budowa grafu deBruijna na podstawie poprawionych odczytów.
4. Usunięcie węzłów stworzonych na podstawie k-merów z listy wrong_kmers, pod warunkiem, że nie są „głową” (węzeł, 
    do którego nie wchodzą żadne krawędzie) albo „ogonem” (węzeł z którego nie wychodzą żadne krawędzie).
    - Końce odczytów mają zazwyczaj niższe pokrycie, więc usuwanie słabych początkowych bądź końcowych k-merów może 
    oznaczać usunięcie poprawnego k-meru pochodzącego z początku bądź końca odczytu. Tego nie chcemy, stąd usuwane są
    tylko słabe k-mery spełniające wymieniony warunek.
5. Usunięcie krawędzi o wsparciu poniżej danego progu, gdzie próg to max(thresh, 1).
    - Thresh to wartość progowa z punktu 2. Można przypuszczać, że słabe krawędzie prowadzą do węzłów utworzonych na 
    podstawie błędnych k-merów. Dzięki takiemu zabiegowi graf jest uproszczony i bardziej pewny.
6. Usuwanie pojedynczych węzłów.
    - Węzeł uznaje się za pojedynczy jeśli nie wchodzi ani nie wychodzi z niego żadna krawędź.
7. Upraszczanie grafu – łączenie węzłów, które są połączone jedynie ze sobą nawzajem.
8. Budowa contigów.
    - Dla każdego węzła-głowy szukane są wszystkiego możliwe contigi, następnie wybierany jest najlepszy z nich (o 
    długości min 300 nukleotydów). W ten sposób w jednej iteracji otrzymujemy po maksymalnie jednym contigu dla 
    każdego węzła początkowego. Z nich wybierany jest jeden najlepszy. Zostaje on zapisany na listę ostatecznych 
    contigów, a następnie wszystkie węzły go tworzące zostają usunięte. Na tak okrojonym zbiorze węzłów 
    przeprowadzana jest kolejna iteracja.
9. Zapis znalezionych contigów do pliku.

<br></br>
<br></br>
<br></br>
#####English

### A tool for assembling reads originating from single-end DNA sequencing.
Instructions for the assignment can be found in assignment2en.pdf.

#### Versions:
- Python 3.6.6

#### How to use?
./assembly input_reads.fasta output_contigs.fasta

#### The process of assembling the reads:
1. Reading reads from fasta file.
2. Error correction: 
    - based on the information about how many and what kind of k-mers can be created from the given reads, mean coverage 
of every k-mer and its standard deviation is received. A difference between mean and standard deviation is taken as the 
threshold. If some k-mer's coverage is below the threshold the other "related" k-mer (Hamming distance equal to 1) is 
searched for. If coverage of related k-mer is above the threshold, the input k-mer is changed into the related one. K-mers 
which coverage is lower than threshold and for which there is no enough good related k-mer are written to a list named 
"wrong_kmers".
3. Construction of DeBruijn graph based on the corrected reads.
4. Deletion of nodes which are based on k-mers from "wrong_kmers" list provided that the node is not a "head" (the node
 without edges coming in) or a "tail" (the node without edges coming out).
    - Ends of the reads usually have the lower coverage, so deletion of "weak" head-like or tail-like k-mers might causes 
    getting rid of the correct k-mers derived from the beginning or the end of the read. We do not want that, so 
    only not head-like or not tail-like nodes based on the weak k-mers are deleted.
5. Removal of the edges which value is below a max(threshold, 1).
    - Threshold is the value from the second point. We can assume that the "weak" edges (which value is low) lead to the
     nodes based on the k-mers with an error. Thanks to such a manipulation the graph is more certain.
6. Deletion of the individual nodes.
    - A node is called individual if no edges come into it or come out of it.
7. Simplifying of the graph - pair of nodes which are connected only to each other is changed into one bigger node.
8. Contigs construction.
     - For each head-node all possible contigs are searched for. Next the best of them (which is longer than 300 
     nucleotides) is chosen. In this way in every iteration for every head one contig is obtained. Then the best of them
      is chosen - it is written to a list of the final contigs. Subsequently all nodes which built the best contig are 
      removed. The next iteration is based on such a truncated set of nodes.
9. Saving found contigs to a file.
