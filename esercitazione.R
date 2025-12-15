#dataset prova


# Carica il pacchetto Biobase, che fornisce le classi base
# per oggetti biologici (es. ExpressionSet) usate in Bioconductor
library(Biobase)



# Carica il pacchetto NMF, che implementa la Non-negative Matrix Factorization
# e funzioni correlate per l’analisi di matrici non negative (es. espressione genica)
library(NMF)


setwd("C:/Users/flavi/OneDrive - Università degli Studi di Bari/_ROOT02_to_update/lavori/conferenze-corsi/2025-12-Bari_workshop_prinn/")
#setwd("C:/Users//Users/nicolettadelbuono/Library/CloudStorage/OneDrive-UniversitàdegliStudidiBari/cartelle_nicoletta_Mac/PresentazioniVarie_&_Convegni/2DaysonCAIMOD-11-12-dic-25/tutorial11-12-25/")

# Legge il file con le informazioni sui campioni (metadati) da un CSV separato da ';'
# 'header = TRUE' indica che la prima riga contiene i nomi delle colonne
Samples_Info <- read.table("Info_Serie_Combinate.csv", header = TRUE, sep = ";")

# Legge la matrice di espressione (o dati) da un file di testo separato da tab ('\t')
# 'header = TRUE' indica che la prima riga contiene i nomi delle colonne (campioni)
# 'row.names = 1' usa la prima colonna come nomi di riga (es. geni)
D <- as.matrix(read.table("Dataset_prova_workshop.txt",
                          header = TRUE,
                          row.names = 1,
                          sep = "\t"))

# Salva i nomi delle righe della matrice D (es. simboli genici) in un vettore
gene_symb <- row.names(D)

# Carica il pacchetto limma, utile per analisi di espressione e diagnostica
library(limma)

# Disegna le densità delle intensità/valori di espressione per tutte le colonne di D
# (senza legenda) per verificare la distribuzione dei dati
plotDensities(D, legend = FALSE)



Wildtype_res<-Samples_Info$WildtypeVSres


# Mostra la dimensione della matrice D: numero di geni (righe) e campioni (colonne)
dim(D)

# Ridisegna le densità, colorando i campioni in base al gruppo Wildtype/resistente
# (funzione 'group' usa il vettore di etichette Wildtype_res)
plotDensities(D, group = Wildtype_res, legend = FALSE)




# PCA
# Trasposta della matrice A: ora le righe sono le osservazioni, le colonne le variabili
X <- t(D)

# Struttura dell'oggetto X (dimensioni, tipo, ecc.)
str(X)

# n = numero di osservazioni (righe), m = numero di variabili (colonne)
n <- nrow(X)
m <- ncol(X)

# Vettore delle medie di ciascuna variabile (colonna)
mean_X <- colMeans(X)

# Costruisce una matrice n x m in cui ogni riga è il vettore delle medie
# (serve per centrare le osservazioni sottraendo la media di ciascuna variabile)
matrix_mean_X <- matrix(mean_X, nrow = n, ncol = m, byrow = TRUE)

# Controlla le prime righe della matrice delle medie
head(matrix_mean_X)

# Centra i dati: sottrae la media di ogni variabile a tutte le osservazioni
X_centered <- X - matrix_mean_X

# Verifica che le medie delle variabili centrate siano ~0 (può esserci solo rumore numerico)
colMeans(X_centered)


# PCA
# Matrice di correlazione delle variabili centrate
# (equivalente a fare PCA sulla matrice standardizzata)
Corr_X <- cor(X_centered)
Corr_X

# Autovalori e autovettori della matrice di correlazione
A <- eigen(Corr_X)
str(A)

# Matrice degli autovettori (ogni colonna è un autovettore = una componente principale)
A_autovect <- A$vectors
str(A_autovect)

# Vettore degli autovalori associati alle componenti principali
lambda <- A$values
lambda

# Percentuale di varianza spiegata da ciascuna componente
# (per la matrice di correlazione la traccia = m)
perc_var_expl <- lambda / m

# Varianza spiegata cumulata (somma progressiva delle percentuali)
var_cum_perc <- cumsum(perc_var_expl)


# Scree plot degli autovalori
plot(lambda,
     type = "b",
     main = "Scree Plot Iris",
     xlab = "Number Component",
     ylab = "Eigenvalue")

# Kaiser's rule:
# Mantieni le componenti con autovalore > 1 (per PCA su matrice di correlazione)
plot(lambda,
     type = "b",
     main = "Scree Plot Iris",
     xlab = "Number Component",
     ylab = "Eigenvalue")
abline(h = 1, lwd = 3, col = "red")

# Dimensione della matrice degli autovettori (m x m)
dim(A_autovect)

# Interpretazione PCA
# Selezione delle prime 2 componenti principali (taglio delle componenti)
A_autovect_red <- A_autovect[, 1:2]
dim(A_autovect_red)

# Proiezione dei dati centrati sul sottospazio delle prime 2 componenti
# (coordinate delle osservazioni nello spazio delle componenti principali)
Y <- X_centered %*% A_autovect_red
dim(Y)


Samples_Info$WildtypeVSre=as.factor(Samples_Info$WildtypeVSre)

dev.new()
plot(Y,main="Bi-plot 2D", xlab = "PC1",ylab="PC2",col = Samples_Info$WildtypeVSre)
abline(v = 0, h=0,col="blue")  


###########################
# PCA con prcomp()
out_pca_prcomp <- prcomp(X_centered, scale = TRUE)

# Riepilogo dell’oggetto PCA: deviazioni standard, proporzione di varianza,
# varianza cumulata per ciascuna componente principale
summary(out_pca_prcomp)
Summary_var <- summary(out_pca_prcomp)

# Loadings (rotazioni): coefficienti che definiscono le componenti principali
# (ogni colonna è una PC, ogni riga è una variabile originale)
head(out_pca_prcomp$rotation)

# Confronto (solo per controllo) con la matrice di autovettori calcolata a mano
head(A_autovect)

# Autovalori della PCA: sdev^2 sono le varianze spiegate da ogni componente
eigvals <- (out_pca_prcomp$sdev)^2
head(eigvals)
# autovalori calcolati prima dalla matrice di correlazione
head(lambda)

# Percentuale di varianza spiegata da ogni componente
perc_pca_prcomp <- eigvals / sum(eigvals)
head(perc_pca_prcomp)

# Controllo (di nuovo) del riepilogo
summary(out_pca_prcomp)

# Scree plot “manuale” delle percentuali di varianza
plot(perc_pca_prcomp, type = "b",
     xlab = "Component",
     ylab = "Proportion of variance",
     main = "Variance explained (prcomp)")

# Varianza spiegata cumulata
cumsum(perc_pca_prcomp)

# Riepilogo (ancora) per vedere varianza cumulata da summary()
summary(out_pca_prcomp)

# Scree plot della varianza cumulata
plot(cumsum(perc_pca_prcomp), type = "b",
     xlab = "Component",
     ylab = "Cumulative proportion",
     main = "Cumulative variance (prcomp)")

# Screeplot “built-in” di R per oggetto prcomp
screeplot(out_pca_prcomp)

# Coordinate delle osservazioni nello spazio delle PC (scores)
X_rec <- out_pca_prcomp$x
dim(X_rec)

# Scatterplot delle prime due componenti principali (PC1 vs PC2)
plot(X_rec[, 1], X_rec[, 2],
     xlab = "PC1",
     ylab = "PC2",
     main = "Scores su PC1–PC2")
abline(v = 0, h = 0, col = "blue")

# Visualizzazione con ggplot2, colorando per tipo di resistenza
library(ggplot2)

# Trasforma le coordinate delle PC in data frame
data_frame_rot <- as.data.frame(out_pca_prcomp$x)

# Aggiunge al data frame la variabile di resistenza dal data frame dei campioni
data_frame_rot$tipo_resistenza <- Samples_Info$tipo_resistenza

head(data_frame_rot)

# Scatterplot PC1–PC2 colorato per Wildtype vs resistente
p <- ggplot(data_frame_rot,
            aes(x = PC1, y = PC2, color = Samples_Info$WildtypeVSres)) +
  geom_point()

# Apre una nuova finestra grafica (utile da Rgui/RStudio su alcune piattaforme)
dev.new()

# Mostra il grafico
p

# Calcola la percentuale di varianza spiegata da ciascuna PC e la esprime in %
percentage <- round(perc_pca_prcomp * 100, 2)

# Costruisce le etichette degli assi: es. "PC1 (45.32 %)"
percentage <- paste(colnames(data_frame_rot),
                    "(",
                    paste(as.character(percentage), "%", ")"))

# Scatter PC1 vs PC2 colorato per tipo di resistenza
p12 <- ggplot(data_frame_rot,
              aes(x = PC1, y = PC2, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[1]) +  # etichetta PC1 con % varianza
  ylab(percentage[2])    # etichetta PC2 con % varianza

# Scatter PC1 vs PC3
p13 <- ggplot(data_frame_rot,
              aes(x = PC1, y = PC3, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[3])

# Scatter PC1 vs PC4
p14 <- ggplot(data_frame_rot,
              aes(x = PC1, y = PC4, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[4])

# Scatter PC1 vs PC5
p15 <- ggplot(data_frame_rot,
              aes(x = PC1, y = PC5, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[5])

# Scatter PC2 vs PC3
p23 <- ggplot(data_frame_rot,
              aes(x = PC2, y = PC3, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[2]) +
  ylab(percentage[3])

# Scatter PC2 vs PC4
p24 <- ggplot(data_frame_rot,
              aes(x = PC2, y = PC4, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[2]) +
  ylab(percentage[4])

# Scatter PC2 vs PC5
p25 <- ggplot(data_frame_rot,
              aes(x = PC2, y = PC5, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[2]) +
  ylab(percentage[5])

# Scatter PC3 vs PC4
p34 <- ggplot(data_frame_rot,
              aes(x = PC3, y = PC4, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[3]) +
  ylab(percentage[4])

# Scatter PC3 vs PC5
p35 <- ggplot(data_frame_rot,
              aes(x = PC3, y = PC5, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[3]) +
  ylab(percentage[5])

# Scatter PC4 vs PC5
p45 <- ggplot(data_frame_rot,
              aes(x = PC4, y = PC5, color = tipo_resistenza)) +
  geom_point() +
  xlab(percentage[4]) +
  ylab(percentage[5])

# Carica il pacchetto grid, che fornisce le funzioni grafiche di basso livello
# su cui si appoggiano sistemi come lattice e ggplot2
library(grid)

# Carica il pacchetto gridExtra, che estende grid con funzioni ad alto livello
# per combinare più grafici (es. ggplot) sulla stessa pagina (grid.arrange, ecc.)
library(gridExtra)

# Carica il pacchetto ggplot2 per la grafica basata su grammatica dei grafici
library(ggplot2)

# Controlla che la variabile 'tipo_resistenza' sia già presente nel data frame
# Se non esiste, la crea copiando i valori dalla colonna corrispondente in Samples_Info
if (!"tipo_resistenza" %in% names(data_frame_rot)) {
  data_frame_rot$tipo_resistenza <- Samples_Info$tipo_resistenza
}

# Rimuove le legende da tutti i plot delle PC (utile se si vogliono combinare più grafici
# in una griglia e mostrare la legenda solo una volta altrove)
p12  <- p12  + theme(legend.position = "none")
p13  <- p13  + theme(legend.position = "none")
p14  <- p14  + theme(legend.position = "none")
p15  <- p15  + theme(legend.position = "none")
p23  <- p23  + theme(legend.position = "none")
p24  <- p24  + theme(legend.position = "none")
p25  <- p25  + theme(legend.position = "none")
p34  <- p34  + theme(legend.position = "none")
p35  <- p35  + theme(legend.position = "none")
p45  <- p45  + theme(legend.position = "none")

# Converte ogni oggetto ggplot in un grob (grafico “pronto” per grid)
# Usa tryCatch per intercettare eventuali errori di conversione
grobs <- lapply(
  list(p12, p13, p14, p15, p23, p24, p25, p34, p35, p45),
  function(p) tryCatch(ggplotGrob(p), error = function(e) e)
)

# Verifica se qualcuno degli elementi è un oggetto di classe "error"
# (cioè la conversione in grob è fallita)
errs <- sapply(grobs, function(x) inherits(x, "error"))

# Se almeno un plot non può essere convertito in grob, interrompe l’esecuzione
# e stampa gli indici dei plot problematici
if (any(errs)) stop("Almeno un plot non può essere convertito in grob. Indici con errore: ",
                    paste(which(errs), collapse = ", "))

# Definisce un grob vuoto (nullGrob) da usare come “segnaposto” per celle vuote in una griglia
empty <- grid::nullGrob()

# Ordine come vuoi (4 righe x 5 colonne):
layout_grobs <- list(
  grobs[[1]], empty,      empty,      empty,      empty,   # p12
  grobs[[2]], grobs[[5]], empty,      empty,      empty,   # p13 p23
  grobs[[3]], grobs[[6]], grobs[[8]], empty,      empty,   # p14 p24 p34
  grobs[[4]], grobs[[7]], grobs[[9]], grobs[[10]], empty    # p15 p25 p35 p45
)

# disegna
grid.newpage()
do.call(grid.arrange, c(layout_grobs, ncol = 5))




######################################
#NMF
######################################


# ora le righe sono le variabili (es. geni) e le colonne sono i campioni,
# come richiesto in genere da molte implementazioni NMF (features x samples)
X_NMF <- D



# (Nota nel codice) la matrice può essere grande; si può mostrare solo una parte in aula
# Stima del rango ottimale (numero di fattori) usando il metodo di Brunet
# Si prova un range di ranghi da 2 a 10, con 2 esecuzioni (nrun = 2) per ciascun rango
# 'seed' serve per rendere riproducibili i risultati
rank_brunet <- nmfEstimateRank(X_NMF, 2:6, method = "brunet", nrun = 50, seed = 123456)

# Grafico delle metriche di qualità (es. consenso) in funzione del rango
# utile per scegliere un numero di fattori ragionevole
plot(rank_brunet)

# Estrae la lista delle matrici di consenso dallo stimatore del rango
cs <- rank_brunet$consensus

# Mappa di consenso per alcuni ranghi (qui ad esempio elementi 1 e 3 della lista)
# con annotazione delle colonne in base a Wildtype vs resistente
consensusmap(cs[c(1, 3)], annCol = Samples_Info$WildtypeVSres)

# Mappa di consenso per altri ranghi (3 e 4) con annotazione per Serie sperimentale
consensusmap(cs[3:4], annCol = Samples_Info$Series)


# Fattorizzazione NMF: scelta del rango
r <- 2

# Esegue una NMF con metodo 'brunet', rango r e una sola esecuzione (nrun = 1)
# .opt = 'vP8' abilita alcune opzioni (es. verbose, parallel, ecc. a seconda della macchina)
out.multirun <- nmf(X_NMF, r, method = "brunet", nrun = 1, .opt = "vP8")

# Valuta il fitting del modello NMF (metriche di errore/adeguatezza)
fit(out.multirun)

# Riepilogo della soluzione NMF (misure di qualità, dimensioni, ecc.)
summary(out.multirun)

# Stampa un breve riassunto dell’oggetto NMF (rango, metodo, ecc.)
out.multirun

# Ricostruzione (approssimata) della matrice originale: W %*% H
out.multirun.cap <- fitted(out.multirun)
dim(out.multirun.cap)
head(out.multirun.cap)

# Matrice dei pesi dei geni (basis): W (features x fattori)
W.tot <- basis(out.multirun)
dim(W.tot)
head(W.tot)

# Visualizza una heatmap delle basi (pattern di espressione dei geni per ciascun fattore)
dev.new()
basismap(out.multirun)

# Matrice dei coefficienti per i campioni: H (fattori x campioni)
H.tot <- coef(out.multirun)
dim(H.tot)
head(H.tot)

# Visualizza una heatmap dei coefficienti, con annotazione per Wildtype vs resistente
# (utile per vedere se i fattori separano i gruppi biologici)
dev.new()
coefmap(out.multirun, annCol = Samples_Info$WildtypeVSres)







#######################################
#clustering
#######################################


# CLUSTERING GERARCHICO

# help su dist: spiega le diverse distanze disponibili (euclidea, manhattan, ecc.)
help(dist)

# Calcola la matrice delle distanze euclidee tra le righe di X
# (qui X ha le osservazioni sulle righe, quindi dist misura distanza tra campioni)
D_euclidea <- dist(X)

# Scalare
# ATTENZIONE: dist() restituisce un oggetto di classe "dist", non una matrice piena;
# il transpose t(D_euclidea) non è davvero necessario/standard.
# Qui si applica una standardizzazione (z-score) alle distanze per avere media ~0 e sd ~1,
# ma spesso per hclust si usa direttamente D_euclidea.
D_euclidea_scaled <- scale(t(D_euclidea))

# Controllo della media e deviazione standard delle distanze scalate
mean(D_euclidea_scaled)
sd(D_euclidea_scaled)

# help su hclust: metodi di collegamento (complete, average, single, ward.D, ecc.)
help(hclust)

# Clustering gerarchico con metodo "complete linkage" sulla distanza euclidea
Clust_G_Compl <- hclust(D_euclidea)

# Struttura dell'oggetto dendrogramma (merge, height, labels, ecc.)
str(Clust_G_Compl)

# Visualizzazione del dendrogramma
dev.new()
plot(Clust_G_Compl)

# Plot con le foglie allineate (hang = -1)
plot(Clust_G_Compl, hang = -1)

# Disegna una linea orizzontale a metà dell'altezza massima del dendrogramma
# (esempio di altezza di taglio per ottenere un certo numero di cluster)
abline(h = max(Clust_G_Compl$height) / 2, col = "red")



# Cambio metodo di linkaggio: average (UPGMA)
# Invece di "complete", si usa la media delle distanze tra gruppi
Clust_G_Av <- hclust(D_euclidea, method = "average")

# Struttura dell’oggetto di clustering gerarchico con linkaggio average
str(Clust_G_Av)

# Nuova finestra grafica e plot del dendrogramma con metodo average
dev.new()
plot(Clust_G_Av)

# Linea di taglio da regolare in base alla scala del dendrogramma
# per suggerire un possibile numero di cluster
abline(h = max(Clust_G_Av$height) / 2, col = "red")

#Tagliare
labels_G_Compl = cutree(Clust_G_Compl,k=2)
table(labels_G_Compl)

#info esterne
table(labels_G_Compl,Samples_Info$WildtypeVSres)


###########
# K-means
##########

# help su kmeans: spiega parametri come 'centers', 'nstart', 'iter.max', algoritmo, ecc.
help("kmeans")

# Esegue il clustering k-means sui campioni (righe di X) chiedendo 2 cluster
# 'centers = 2' specifica il numero di cluster; conviene spesso usare nstart > 1
Clust_K <- kmeans(X, centers = 2)

# Struttura dell’oggetto kmeans: cluster assegnati, centri, withinss, tot.withinss, ecc.
str(Clust_K)

#label,centroidi,cardinalit?
labels_K = Clust_K$cluster
Clust_K$centers
Clust_K$size

#label note
table(labels_K,Samples_Info$WildtypeVSres)



# Numero ottimale di cluster

# Metodo WSS (Within-Cluster Sum of Squares) “a mano”
Kmax <- 9                 # numero massimo di cluster da provare
wit_ss <- NULL            # vettore per salvare il WSS totale per ogni k

for (i in seq(1, Kmax)) {
  # Esegue k-means con i cluster e salva la somma delle varianze intra-cluster
  W_aux <- kmeans(X, centers = i)$tot.withinss
  wit_ss[i] <- W_aux
}

wit_ss
length(wit_ss)

# Grafico WSS vs numero di cluster (metodo del “gomito”)
plot(1:Kmax, wit_ss,
     xlab = "Numero Cluster",
     ylab = "Tot within-cluster SS",
     main = "WSS",
     type = "b")   # linee + punti

# Uso di factoextra per la scelta del numero di cluster

library(factoextra)
help(fviz_nbclust)

# Metodo WSS/elbow automatizzato
fviz_nbclust(X, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +      # esempio di gomito a k = 4
  labs(subtitle = "WSS (metodo del gomito)")

# Metodo della silhouette media
fviz_nbclust(X, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette")

###########
#misure di Validazione del Clustering
###########
#Clustering Gerarchico 
#Clust_G_Compl

#Silhouette
myclusters = cutree(Clust_G_Compl,2)
table(myclusters,Samples_Info$WildtypeVSres)

library(cluster)
Sil_coef = silhouette(myclusters,D_euclidea)
str(Sil_coef)
plot(Sil_coef,main = "Gerarchico Complete Euclidea",border=NA,col=1:3)

#proviamo a cambiare cut
myclusters = cutree(Clust_G_Compl,2)
# table(myclusters,iris_Species)
# 
Sil_coef = silhouette(myclusters,D_euclidea)
str(Sil_coef)
plot(Sil_coef,main = "Gerarchico Complete Euclidea",border=NA,col=1:4)

#Estrapolare info da Sil Coef
summary_sil = summary(Sil_coef)
str(summary_sil)

summary_sil$clus.avg.widths
summary_sil$avg.width
summary_sil$clus.sizes


library(factoextra)
fviz_silhouette(Sil_coef)
