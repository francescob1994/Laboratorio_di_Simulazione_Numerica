- main file = Monte_Carlo_NVT.cpp

per il calcolo delle correlazioni:
- repeat.sh esegue 60 simulazioni inizializzando l'oggetto random in modo diverso ad ogni simulazione. La simulazione n-esima
  stampa nel file u_p.n valori istantanei di energia e pressione. In seguito viene eseguito lo script correlazioni.py con input 'n'.
- correlazioni.py stampa sul file cov.n nella prima colonna i prodotti delle energie istantanee con la prima energia e nella seconda 
  gli stessi prodotti per la pressione
- mean.sh genera il file up.mean che contiene le medie delle energie e pressioni istantanee e dei prodotti calcolati con correlazioni.py 
  mediando sulle 60 simulazioni
- fit.py calcola l'andamento delle correlazioni di energia e pressione e calcola la lunghezza di correlazione

- differentL.sh lancia 30 simulazioni (inizializzando l'oggetto random in modo diverso ad ogni simulazione) e viene utilizzato 
  nel calcolo dell'andamento dell'errore con la dimensione dei blocchi
- gr.py disegna le funzioni g(r)  
