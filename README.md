# Tutorial: Fattorizzazioni low-rank su matrici di dati omici e clustering

## Descrizione

Questo repository contiene il materiale del tutorial **"Fattorizzazioni low-rank su matrici di dati omici e clustering"**, pensato come introduzione pratica all’uso di tecniche di fattorizzazione di matrici e clustering in **R** applicate a dati omici e clinici.

Il tutorial utilizza come caso studio un dataset tratto dalla letteratura e guida passo passo nell’applicazione di:

* **Principal Component Analysis (PCA)**
* **Nonnegative Matrix Factorization (NMF)**
* **Clustering gerarchico**
* **k-means**

con particolare attenzione all’interpretazione e alla visualizzazione dei risultati.

---

## Evento

**Two days on Computational Approaches for the Integration of Multi-Omics Data**
Aula VII, Department of Mathematics – University of Bari Aldo Moro
**11–12 dicembre 2025**

Sito web dell’evento:
[https://sites.google.com/view/caimod-projectprin-2022pnrr/two-days-on-caimod](https://sites.google.com/view/caimod-projectprin-2022pnrr/two-days-on-caimod)

---

## Relatori

* **Nicoletta Del Buono**
* **Flavia Esposito**

Dipartimento di Matematica, Università degli Studi di Bari Aldo Moro

---

## Obiettivi di apprendimento

Al termine del tutorial, i partecipanti saranno in grado di:

1. Comprendere i concetti fondamentali delle **fattorizzazioni di matrici** e del **clustering**.
2. Applicare i principali **metodi numerici** associati a tali tecniche.
3. Utilizzare e interpretare **metodi di visualizzazione** per spiegare e comunicare i risultati ottenuti.

---

%## Basi teoriche

%Il tutorial mette in pratica gli argomenti presentati nella parte teorica.
%Le slide della parte teorica sono disponibili nel file:


---

## Struttura del repository

Una possibile organizzazione dei file è la seguente:

```
Tutorial/
│── README.md
│── esercitazione.R
│── data/
│   └── dataset_omico.csv
```

*(La struttura può variare in base ai file effettivamente distribuiti durante il tutorial.)*

---

## Istruzioni pratiche

Per svolgere il tutorial è necessario installare:

* **R**, disponibile per il proprio sistema operativo dal sito CRAN (IT Mirror consigliato)
* **RStudio**, l’IDE distribuito da **Posit**

### Download dei file

1. Scaricare tutti i file presenti in questo repository.
2. Creare una cartella (ad esempio `Tutorial`).
3. Posizionare **tutti i file del tutorial nella stessa working directory di R**.

> **Nota:** è fortemente consigliato lavorare all’interno di una singola cartella per evitare problemi con i percorsi dei file.

---

## Ringraziamenti

Questo tutorial fa parte delle attività di disseminazione del progetto:

**Computational approaches for the integration of multi-omics data**
Progetto **P2022BLN38**, finanziato dall’**Unione Europea – Next Generation EU** nell’ambito del programma **PRIN 2022 PNRR**
(D.D. 1409 del 14-09-2022, Ministero dell’Università e della Ricerca)

**CUP:** B53D23027810001

This tutorial is part of the dissemination activities supported by the P2022BLN38 project Computational approaches for the integration of multi-omics data – funded by European Union – Next Generation EU within the PRIN 2022 PNRR program (D.D. 1409 del 14-09-2022 Ministero dell’Università e della Ricerca) CUP B53D23027810001.
---

## Contatti

Per domande o segnalazioni, aprire una *issue* su GitHub oppure contattare gli autori del tutorial.
