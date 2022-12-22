# Problem n teles
## Opis
Mnoge probleme v fiziki in kemiji lahko rešimo s proučevanjem interakcij med vsemi
telesi v sistemu. Na področju molekulske dinamike je na primer zanimivo obnašanje
molekul v Lennard-Jonesovem potencialu, na področju astronomije pa vpliv gravitacijskih
sil na gibanje teles. Ker računamo interakcije vseh teles na izbrano telo, je računska
kompleksnost osnovnega sekvenčnega algoritma O(n2) za vsak časovni korak, kjer je n
število vseh teles v sistemu. Kljub temu, da obstaja mnogo izboljšav osnovnega algoritma,
se bomo v nadaljevanju omejili na vzporejanje osnovnega algoritma.

## Naloga
Pripravite program za sekvenčno simulacijo obnašanja sistema več nebesnih teles. Gra-
vitacijsko kontantko κ in mase teles nastavite tako, da bo potek simulacije zanimiv.
Program prilagodite za učinkovito računanje na sistemih s skupnim in sistemih s poraz-
deljenim pomnilnikom. Rešitev mora delovati za poljubno število teles in za poljubno
število procesorjev. Nekaj idej:

- Primerjajte vzporejanje s knjižnicama Pthread in OpenMP.
- Problem lahko razdelite med procesorje tako, da vsak računa skupno silo na izbrano število teles. Ker pa velja Fij = −Fij , lahko število izračunov prepolovite.Primerjajte obe rešitvi.
- Poskrbite za učinkovito shranjevanje položaja teles po vsakem k-tem koraku simulacije, mogoče MPI-IO in izdelajte film.
- Pri MPI preverite ali je podatke o telesih bolje organizirati v polje struktur (podatkov o telesih) ali v strukturo z več polji (eno za položaj, drugo za hitrost vseh teles).
- Raziščite uporabo podatkovnih tipov MPI.
- Za prenašanje podatkov med vozlišči uporabite knjižnico MPI, za računanje na istem vozlišču pa knjižnico OpenMP.
- Uporabite asinhrono prenašanje podatkov.
