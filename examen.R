

install.packages("ChainLadder")
library(ChainLadder)

#DATOS DESACOMULADOS
c0<-c(4946974,3721237,894717,207760,206704,62124,64813,14840,11130,14813)
c1<-c(6346746,3246406,723222,141797,67824,36603,42742,11186,11646)
c2<-c(6269090,2976223,847043,262768,152703,65444,53545,8924)
c3<-c(6863016,2683224,722532,190653,132976,88340,43329)
c4<-c(6778886,2746229,653894,273395,230288,105224)
c5<-c(6184793,2828338,572765,244899,104957)
c6<-c(6600184,2893207,563114,225517)
c7<-c(6288066,2440103,528043)
c8<-c(5290793,2357936)
c9<-c(5675568)

vpf<-ibnrchl(c(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9))
sum(vpf)

#CONSIDERANDO ETTI
ETI.V<-c(0.00024,0.00031,0.00052,0.00120,0.00175,0.00749,0.00963,0.00985,0.00987)
i.renta<-numeric(length(ETI.V))
for (i in 1:length(ETI.V)) {i.renta[i]<-(1+ETI.V[i])^(-i)}
prov<-sum(vpf*i.renta);prov

#CON FLUJOS
tf<-numeric(length(ETI.V))
tf[1]<-ETI.V[1]
for (i in 2:length(ETI.V)) {tf[i]<-(((1+ETI.V[i])^(i))/((1+ETI.V[i-1])^(i-1)))-1}

f.fluxe<-numeric(length(ETI.V))
for (i in 1:length(ETI.V)) { f.fluxe[i]<-(1-(1+tf[i])^(-1))/(log(1+tf[i]))}

i.fluxe<-rep(1,length(ETI.V))
for (i in 2:length(ETI.V)) {i.fluxe[i]<-(1+ETI.V[i-1])^(-i+1)}

provi<-sum(vpf*i.fluxe*f.fluxe);provi

#MODELO DE MACK
C0<-cumsum(c0)
C1<-c(cumsum(c1),NA)
C2<-c(cumsum(c2),NA,NA)
C3<-c(cumsum(c3),NA,NA,NA)
C4<-c(cumsum(c4),NA,NA,NA,NA)
C5<-c(cumsum(c5),NA,NA,NA,NA,NA)
C6<-c(cumsum(c6),NA,NA,NA,NA,NA,NA)
C7<-c(cumsum(c7),NA,NA,NA,NA,NA,NA,NA)
C8<-c(cumsum(c8),NA,NA,NA,NA,NA,NA,NA,NA)
C9<-c(cumsum(c9),NA,NA,NA,NA,NA,NA,NA,NA,NA)

C<-matrix(c(C0,C1,C2,C3,C4,C5,C6,C7,C8,C9),ncol = 10)
C<-t(C); C<-as.triangle(C)

mhc<-MackChainLadder(C);mhc

#MACK
desac<-cum2incr(mhc$FullTriangle)
vpf<-numeric(dim(C)[1]-1)

for (k in 1:(dim(C)[1] - 1)) 
  { futuro<-(row(desac)+col(desac)-1)==(dim(C)[1]+k)
  vpf[k]<-sum(desac[futuro])
}

proveti<-sum(vpf*i.renta);proveti
proflux<-sum(vpf*i.fluxe*f.fluxe);proflux

#POISSON
poi<-glmReserve(C,var.power = 1,link.power = 0,mse.method = "formula");poi

desac<-cum2incr(poi$FullTriangle)
vpf<-numeric(dim(C)[1]-1)

for (k in 1:(dim(C)[1] - 1)) 
{ futuro<-(row(desac)+col(desac)-1)==(dim(C)[1]+k)
vpf[k]<-sum(desac[futuro])
}

proveti<-sum(vpf*i.renta);proveti
proflux<-sum(vpf*i.fluxe*f.fluxe);proflux


#GAMMA
gamma<-glmReserve(C,var.power = 2, link.power = 0, mse.method = "formula");gamma

desac<-cum2incr(gamma$FullTriangle)
vpf<-numeric(dim(C)[1]-1)

for (k in 1:(dim(C)[1] - 1)) 
{ futuro<-(row(desac)+col(desac)-1)==(dim(C)[1]+k)
vpf[k]<-sum(desac[futuro])
}

proveti<-sum(vpf*i.renta);proveti
proflux<-sum(vpf*i.fluxe*f.fluxe);proflux





# ------------------------------------------------------------------------------
# Supuesto 1
# ------------------------------------------------------------------------------

library(ChainLadder)

load("provisio.RData")

c0 <- c(88,43.6,51,54.15,15.6)
c1 <- c(93.2,45,64.2,54.8)
c2 <- c(109,69.2,57.4)
c3 <- c(122.4,63.4)
c4 <- 136.8

ibnrchl(c(c0,c1,c2,c3,c4))

ibnrchlvar1(c(c0,c1,c2,c3,c4))

ibnrchlvar2(c(c0,c1,c2,c3,c4))

ibnrvylder(c(c0,c1,c2,c3,c4))

ni <- c(630,750,800,805,935)

ibnrarit(c(c0,c1,c2,c3,c4),ni)

ibnrgeom(c(c0,c1,c2,c3,c4),ni)

ibnrreg(c(c0,c1,c2,c3,c4),ni)

C0 <- cumsum(c0)
C1 <- c(cumsum(c1),NA)
C2 <- c(cumsum(c2),NA,NA)
C3 <- c(cumsum(c3),NA,NA,NA)
C4 <- c(c4,NA,NA,NA,NA)

C <- matrix(c(C0,C1,C2,C3,C4),ncol=5)
C <- t(C); C
C <- as.triangle(C); C
mch <- MackChainLadder(C); mch
names(mch)
mch$FullTriangle
mch$f

glmformula<-glmReserve(C,var.power = 1, link.power = 0, mse.method = "formula"); glmformula
names(glmformula)
glmformula$FullTriangle
glmformula$model



# ------------------------------------------------------------------------------
# Supuesto 2
# ------------------------------------------------------------------------------

# a) -------------------------------------------------
# a.1) -----------------------------------------------

load("provisio.RData")

c0<-c(3487.2,3003.1,1442.4,803,330)
c1<-c(2972.2,3000.8,1155.3,955.6)
c2<-c(3312.4,3059.6,1431.9)
c3<-c(3894.7,3991.9)
c4<-4483.6
vpf<-ibnrchl(c(c0,c1,c2,c3,c4))

# a.2) -----------------------------------------------

vpf
sum(vpf)
ETI.V<- c(0.00036,0.00015,0.00067,0.00155)
i.renta<-numeric(length(ETI.V))
for (i in 1:length(ETI.V)) {i.renta[i]<- (1+ETI.V[i])^(-i)}; i.renta
prov.renta<-sum(vpf*i.renta); prov.renta

# a.3) -----------------------------------------------

tf<-numeric(length(ETI.V))
tf[1] <- ETI.V[1]
for (i in 2:length(ETI.V)) {tf[i]<- ((1+ETI.V[i])^(i))/((1+ETI.V[i-1])^(i-1))-1}; tf

f.fluxe<-numeric(length(ETI.V))
for (i in 1:length(ETI.V)) {f.fluxe[i]<- (1-((1+tf[i])^(-1)))/log(1+tf[i])}
f.fluxe

i.fluxe<-rep(1,length(ETI.V))
for (i in 2:length(ETI.V)) {i.fluxe[i]<- (1+ETI.V[i-1])^(-i+1)} ; i.fluxe

prov.fluxe<-sum(vpf*i.fluxe*f.fluxe);prov.fluxe

# b) -------------------------------------------------

vpf1<-ibnrchlvar1(c(c0,c1,c2,c3,c4))

vpf1
sum(vpf1)
prov.renta1<-sum(vpf1*i.renta);prov.renta1
prov.fluxe1<-sum(vpf1*i.fluxe*f.fluxe);prov.fluxe1

vpf2<-ibnrchlvar2(c(c0,c1,c2,c3,c4))

vpf2
sum(vpf2)
prov.renta2<-sum(vpf2*i.renta);prov.renta2
prov.fluxe2<-sum(vpf2*i.fluxe*f.fluxe);prov.fluxe2

# ------------------------------------------------------------------------------
# Supuesto 3
# ------------------------------------------------------------------------------

# a) -------------------------------------------------

load("provisio.RData")

c0<-c(288,132,240,40,50)
c1<-c(398,102,198,102)
c2<-c(530,170,90)
c3<-c(610,190)
c4<-715

vpf<-ibnrchl(c(c0,c1,c2,c3,c4))

# Asumiendo que se pagan completamente durante
# el a?o de ocurrencia y los 4 siguientes

vpf
sum(vpf)
prov.renta<-sum(vpf*i.renta);prov.renta

# La idea si no asumimos ETI.V 

importe.total<-750+857.1429+934.9358+1255.3475+1486.9440;importe.total
importe.pagado<-750+800+790+800+715; importe.pagado
importe.pagado<-sum(c0,c1,c2,c3,c4); importe.pagado
total.pagos<-importe.total-importe.pagado;total.pagos

# Asumiendo que se paga el 85% del importe de los siniestros
# durante el a?o de ocurrencia y los 4 siguientes, es decir,
# que el tri?ngulo de datos supone ?nicamente el 85% de todos 
# los pagos

importe.total85<-importe.total ; importe.total85
importe.pagado85<-importe.pagado;importe.pagado85
vpf85<-vpf ;vpf85

importe.total.nuevo<-importe.total85/(85/100); importe.total.nuevo
total.pagos.nuevos<-importe.total.nuevo-importe.pagado85; total.pagos.nuevos
pagos.finales <-c(0, 0, 0, total.pagos.nuevos-sum(vpf85)); pagos.finales 

pagos.pendientes<-vpf85+pagos.finales; pagos.pendientes

prov.renta.nueva<-sum((pagos.pendientes)*i.renta); prov.renta.nueva

# b) -------------------------------------------------

vpf<-ibnrvylder(c(c0,c1,c2,c3,c4))

vpf
sum(vpf)
prov.renta<-sum(vpf*i.renta);prov.renta

# c) -------------------------------------------------

ni<-c(651,752,825,900,1135)
vpf<-ibnrarit(c(c0,c1,c2,c3,c4),ni)
vpf
sum(vpf)
prov.renta<-sum(vpf*i.renta);prov.renta


# ------------------------------------------------------------------------------
# Supuesto 4
# ------------------------------------------------------------------------------

load("provisio.RData")

c0 <- c(1001,854,568,565,347,148)
c1 <- c(1113,990,671,648,422) 
c2 <- c(1265,1168,800,744) 
c3 <- c(1490,1383,1007)  
c4 <- c(1725,2536) 
c5 <- c(1889)

ibnrchl(c(c0,c1,c2,c3,c4,c5))

# ------------------------------------------------------------------------------
# Modelo de Mack
# ------------------------------------------------------------------------------

## Instalamos el paquete "ChainLadder" y lo cargamos

library(ChainLadder)

C0 <- cumsum(c0)
C1 <- c(cumsum(c1),NA)
C2 <- c(cumsum(c2),NA,NA)
C3 <- c(cumsum(c3),NA,NA,NA)
C4 <- c(cumsum(c4),NA,NA,NA,NA)
C5 <- c(c5,NA,NA,NA,NA,NA)

C <- matrix(c(C0,C1,C2,C3,C4,C5),ncol=6)
C <- t(C); C
C <- as.triangle(C); C
mch <- MackChainLadder(C); mch
names(mch)
mch$Triangle
mch$FullTriangle
mch$f

# En FullTriangle est? el acumulado, lo desacumulamos para 
# obtener el vector de pagos futuros

a<-matrix(c(rep(0,dim(C)[1]),mch$FullTriangle),nrow=dim(C)[1],ncol=dim(C)[1]);a
noncumFullTriangle<-mch$FullTriangle-a; noncumFullTriangle

# Obtenemos el vector de pagos futuros

vpf <- rep(0, dim(C)[1] - 1)
for (k in 1:dim(C)[1] - 1) {
        future <- row(noncumFullTriangle) + col(noncumFullTriangle) - 1 == dim(C)[1] + k
        vpf[k] <- sum(noncumFullTriangle[future])
    }
vpf
ETI.V<- c(0.00036,0.00015,0.00067,0.00155,0.00267)
i.renta<-numeric(length(ETI.V))
for (i in 1:length(ETI.V)) {i.renta[i]<- (1+ETI.V[i])^(-i)}; i.renta
prov.renta<-sum(vpf*i.renta); prov.renta

# ------------------------------------------------------------------------------
# GLM Poisson sobredisperso (f?rmula)
# ------------------------------------------------------------------------------

# glmReserve(triangle, var.power = 1, link.power = 0, cum = TRUE, 
# mse.method = c("formula", "bootstrap"), nsim = 1000, nb = FALSE, ...)

glmformula<-glmReserve(C,mse.method = "formula"); glmformula
names(glmformula)
glmformula$FullTriangle

## Observamos que aqu? nos presenta en el FullTriangle el resultado acumulado,
## lo mismo que ocurre en el modelo de Mack

# Una vez convertido C en un objeto de tipo triángulo podemos utilizar 
# las funciones de conversión cum2incr y incr2cum. En Mack hay que incluir
# el triángulo acumulado y en glm podemos elegir, aunque por defecto cum = TRUE

Cincr <- cum2incr(C)
Ccum <- incr2cum(Cincr)

# Lo mismo para el resultado glmformula$FullTriangle:

noncumFullTriangle <- cum2incr(glmformula$FullTriangle); noncumFullTriangle
vpf <- rep(0, dim(C)[1] - 1)
for (k in 1:dim(C)[1] - 1) {
        future <- row(noncumFullTriangle) + col(noncumFullTriangle) - 1 == dim(C)[1] + k
        vpf[k] <- sum(noncumFullTriangle[future])
    }
vpf
prov.renta<-sum(vpf*i.renta);prov.renta

glmformula$model

## Si queremos predecir sin que sean los datos acumulados podemos:
 
# A?o de origen 1 a?o de desarrollo 5
exp(6.78751 )*exp(0.14204)*exp(-1.79030)

# ------------------------------------------------------------------------------
# GLM: Gamma, Logar?tmico
# ------------------------------------------------------------------------------

gammalog <- glmReserve(C,mse.method="formula",var.power=2, link.power = 0)
summary(gammalog)
names(gammalog)
noncumFullTriangle <- cum2incr(gammalog$FullTriangle); noncumFullTriangle
vpf <- rep(0, dim(C)[1] - 1)
for (k in 1:dim(C)[1] - 1) {
        future <- row(noncumFullTriangle) + col(noncumFullTriangle) - 1 == dim(C)[1] + k
        vpf[k] <- sum(noncumFullTriangle[future])
    }
vpf
prov.renta<-sum(vpf*i.renta);prov.renta

# ------------------------------------------------------------------------------
# GLM: Inversa Guassiana, Logar?tmico
# ------------------------------------------------------------------------------

invgausslog <- glmReserve(C,mse.method="formula",var.power=3,link.power=0)
summary(invgausslog)
noncumFullTriangle <- cum2incr(invgausslog$FullTriangle); noncumFullTriangle
vpf <- rep(0, dim(C)[1] - 1)
for (k in 1:dim(C)[1] - 1) {
        future <- row(noncumFullTriangle) + col(noncumFullTriangle) - 1 == dim(C)[1] + k
        vpf[k] <- sum(noncumFullTriangle[future])
    }
vpf
prov.renta<-sum(vpf*i.renta);prov.renta

# ------------------------------------------------------------------------------
# GLM Poisson sobredisperso (bootstrap)
# ------------------------------------------------------------------------------

set.seed(1111)
glmbootstrap1000<-glmReserve(C,mse.method = "bootstrap"); glmbootstrap1000
names(glmbootstrap1000)
glmbootstrap1000$summary
noncumFullTriangle <- cum2incr(glmbootstrap1000$FullTriangle); noncumFullTriangle
vpf <- rep(0, dim(C)[1] - 1)
for (k in 1:dim(C)[1] - 1) {
         future <- row(noncumFullTriangle) + col(noncumFullTriangle) - 1 == dim(C)[1] + k
         vpf[k] <- sum(noncumFullTriangle[future])
     }
vpf
prov.renta<-sum(vpf*i.renta);prov.renta

## glmbootstrap1000$sims.reserve.pred matriz 5 x 1000 con las distribuciones 
## predictivas de las reservas por a?os de origen:

## A?o de origen 1: glmbootstrap1000$sims.reserve.pred[,1]
## A?o de origen 2: glmbootstrap1000$sims.reserve.pred[,2]
## A?o de origen 3: glmbootstrap1000$sims.reserve.pred[,3]
## A?o de origen 4: glmbootstrap1000$sims.reserve.pred[,4]
## A?o de origen 5: glmbootstrap1000$sims.reserve.pred[,5]

mean(glmbootstrap1000$sims.reserve.pred[,1])
mean(glmbootstrap1000$sims.reserve.pred[,2])
mean(glmbootstrap1000$sims.reserve.pred[,3])
mean(glmbootstrap1000$sims.reserve.pred[,4])
mean(glmbootstrap1000$sims.reserve.pred[,5])
 
## Para la distribuci?n predictiva de la reserva total sumamos los valores:

Total<-glmbootstrap1000$sims.reserve.pred[,1]+
glmbootstrap1000$sims.reserve.pred[,2] + 
glmbootstrap1000$sims.reserve.pred[,3] +
glmbootstrap1000$sims.reserve.pred[,4] +
glmbootstrap1000$sims.reserve.pred[,5]

mean(Total)
sd(Total)
cv<-(sd(Total)/mean(Total))*100 ; cv  ## c.v. in %
VaR<-quantile(Total, c(0.5, 0.75, 0.9, 0.95, 0.995)) ; VaR  ## VaR(alfa)
pp<-(Total-mean(Total))/sd(Total)
sum(pp^3)/(1000-1)  # estimates the skewness
sum(pp^4)/(1000-1) -3 # estimates the kurtosis
hist(Total,nclass=20)

# ------------------------------------------------------------------------------
# Supuesto 5
# ------------------------------------------------------------------------------

# Datos acumulados:
C0 <- c(1001,1855,2423,2988,3335,3483)
C1 <- c(1113,2103,2774,3422,3844)
C2 <- c(1265,2433,3233,3977)
C3 <- c(1490,2873,3880)
C4 <- c(1725,4261)
C5 <- 1889
 
fd <- c(2.051107,1.328800,1.232147,1.119969,1.044378)
alphaext <- c(3520,3980,4620,5660,6210,6330)
gammaext <- c(0.28,0.53,0.71,0.86,0.95,1)
gammaest <- rep(0,6)
gammaest[1] <- 1/prod(fd)
gammaest[2] <- 1/prod(fd[2:5])
gammaest[3] <- 1/prod(fd[3:5])
gammaest[4] <- 1/prod(fd[4:5])
gammaest[5] <- 1/prod(fd[5:5])
gammaest[6] <- 1
gammaest

# Caso V11

C15 <- C1[5]+(gammaext[6]-gammaext[5])*alphaext[2];C15
C24 <- C2[4]+(gammaext[5]-gammaext[4])*alphaext[3];C24
C25 <- C2[4]+(gammaext[6]-gammaext[4])*alphaext[3];C25
C33 <- C3[3]+(gammaext[4]-gammaext[3])*alphaext[4];C33
C34 <- C3[3]+(gammaext[5]-gammaext[3])*alphaext[4];C34
C35 <- C3[3]+(gammaext[6]-gammaext[3])*alphaext[4];C35
C42 <- C4[2]+(gammaext[3]-gammaext[2])*alphaext[5]; C42
C43 <- C4[2]+(gammaext[4]-gammaext[2])*alphaext[5];C43
C44 <- C4[2]+(gammaext[5]-gammaext[2])*alphaext[5];C44
C45 <- C4[2]+(gammaext[6]-gammaext[2])*alphaext[5];C45
C51 <- C5[1]+(gammaext[2]-gammaext[1])*alphaext[6];C51
C52 <- C5[1]+(gammaext[3]-gammaext[1])*alphaext[6];C52
C53 <- C5[1]+(gammaext[4]-gammaext[1])*alphaext[6];C53
C54 <- C5[1]+(gammaext[5]-gammaext[1])*alphaext[6];C54
C55 <- C5[1]+(gammaext[6]-gammaext[1])*alphaext[6];C55
Provextext <- C15+C25+C35+C45+C55-C1[5]-C2[4]-C3[3]-C4[2]-C5; Provextext

# Caso V13

C15bis <- C1[5]+(gammaest[6]-gammaest[5])*alphaext[2]; C15bis
C24bis <- C2[4]+(gammaest[5]-gammaest[4])*alphaext[3];C24bis
C25bis <- C2[4]+(gammaest[6]-gammaest[4])*alphaext[3];C25bis
C33bis <- C3[3]+(gammaest[4]-gammaest[3])*alphaext[4];C33bis
C34bis <- C3[3]+(gammaest[5]-gammaest[3])*alphaext[4];C34bis
C35bis <- C3[3]+(gammaest[6]-gammaest[3])*alphaext[4];C35bis
C42bis <- C4[2]+(gammaest[3]-gammaest[2])*alphaext[5];C42bis
C43bis <- C4[2]+(gammaest[4]-gammaest[2])*alphaext[5];C43bis
C44bis <- C4[2]+(gammaest[5]-gammaest[2])*alphaext[5];C44bis
C45bis <- C4[2]+(gammaest[6]-gammaest[2])*alphaext[5];C45bis
C51bis <- C5[1]+(gammaest[2]-gammaest[1])*alphaext[6];C51bis
C52bis <- C5[1]+(gammaest[3]-gammaest[1])*alphaext[6];C52bis
C53bis <- C5[1]+(gammaest[4]-gammaest[1])*alphaext[6];C53bis
C54bis <- C5[1]+(gammaest[5]-gammaest[1])*alphaext[6];C54bis
C55bis <- C5[1]+(gammaest[6]-gammaest[1])*alphaext[6];C55bis
Provestext <- C15bis+C25bis+C35bis+C45bis+C55bis-C1[5]-C2[4]-C3[3]-C4[2]-C5; Provestext

# Caso V11

F0 <- c(1001,1855,2423,2988,3335,3483)
F1 <- c(1113,2103,2774,3422,3844,C15)
F2 <- c(1265,2433,3233,3977,C24,C25)
F3 <- c(1490,2873,3880,C33,C34,C35)
F4 <- c(1725,4261,C42,C43,C44,C45)
F5 <- c(1889,C51,C52,C53,C54,C55)

M<-rbind(F0,F1,F2,F3,F4,F5)
M<-as.matrix(M)
Mt<-t(M)
Mtdiff<-diff(Mt)
Mnoacum<-t(Mtdiff)
Mnoacum<-matrix(c(M[,1],Mnoacum),nrow=6,ncol=6);Mnoacum
vpf <- rep(0, dim(Mnoacum)[1] - 1)
for (k in 1:dim(Mnoacum)[1] - 1) {
        future <- row(Mnoacum) + col(Mnoacum) - 1 == dim(Mnoacum)[1] + k
        vpf[k] <- sum(Mnoacum[future])
    }
vpf
prov.renta<-sum(vpf*i.renta);prov.renta

# Caso V13

F0 <- c(1001,1855,2423,2988,3335,3483)
F1 <- c(1113,2103,2774,3422,3844,C15bis)
F2 <- c(1265,2433,3233,3977,C24bis,C25bis)
F3 <- c(1490,2873,3880,C33bis,C34bis,C35bis)
F4 <- c(1725,4261,C42bis,C43bis,C44bis,C45bis)
F5 <- c(1889,C51bis,C52bis,C53bis,C54bis,C55bis)
 
M<-rbind(F0,F1,F2,F3,F4,F5)
M<-as.matrix(M)
Mt<-t(M)
Mtdiff<-diff(Mt)
Mnoacum<-t(Mtdiff)
Mnoacum<-matrix(c(M[,1],Mnoacum),nrow=6,ncol=6);Mnoacum
vpf <- rep(0, dim(Mnoacum)[1] - 1)
for (k in 1:dim(Mnoacum)[1] - 1) {
        future <- row(Mnoacum) + col(Mnoacum) - 1 == dim(Mnoacum)[1] + k
        vpf[k] <- sum(Mnoacum[future])
    }
vpf
prov.renta<-sum(vpf*i.renta);prov.renta

