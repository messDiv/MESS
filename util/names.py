#!/usr/bin/env python2.7
import random

## Toy function for generating random genus/species names

namelist = """Elateridae  Hypnoidus   abbreviatus (Say)   
   Psyllidae   Acizzia uncatoides  (Ferris & Klyver)   
 Sphingidae  Eumorpha    achemon (Drury) 
   Coccidae    Kilifia acuminata   (Signoret)  
  Cryptophagidae  Cryptophagus    acutangulus Gyllenhal   
 Vespidae    Dolichovespula  arenaria    (Fabricius) 
 Yponomeutidae   Atteva  punctella   (Cramer)    
  Curculionidae   Alniphagus  aspericollis    (LeConte)   
  Chrysomelidae   Altica  ambiens (LeConte)   
   Cercopidae  Clastoptera obtusa  (Say)   
 Pieridae    Colias  eurytheme   Boisduval   
 Cecidomyiidae   Asphondylia websteri    Felt    
 Megachilidae    Megachile   rotundata   (Fabricius) 
 Gelechiidae Dichomeris  acuminata   (Staudinger)    
 Noctuidae   Autographa  californica (Speyer)    
   Miridae Adelphocoris    lineolatus  (Goeze) 
 Eurytomidae Bruchophagus    roddi   (Gussakovsky)   
  Curculionidae   Otiorhynchus    ligustici   (Linnaeus)  
 Crambidae   Loxostege   cereralis   (Zeller)    
  Curculionidae   Hypera  postica (Gyllenhal) 
 Halictidae  Nomia   melanderi   Cockerell   
 Formicidae  Formica exsectoides Forel   
 Pyralidae   Cadra   cautella    (Walker)    
  Chrysomelidae   Gonioctena  americana   (Schaeffer) 
  Tenebrionidae   Tribolium   audax   Halstead    
 Apidae  Bombus  pensylvanicus   (De Geer)   
  Silphidae   Nicrophorus americanus  Olivier 
Blattidae   Periplaneta americana   (Linnaeus)  
 Noctuidae   Acronicta   americana   (Harris)    
 Ixodidae    Dermacentor variabilis  (Say)   
 Sesiidae    Sesia   tibialis    (Harris)    
 Pyroglyphidae   Dermatophagoides    farinae Hughes  
   Cixiidae    Myndus  crudus  Van Duzee   
 Pyralidae   Euzophera   semifuneralis   (Walker)    
 Agromyzidae Liriomyza   trifolii    (Burgess)   
  Ptinidae    Mezium  americanum  (Castelnau) 
 Gelechiidae Sitotroga   cerealella  (Olivier)   
   Cicadellidae    Acinopterus angulatus   Lawson  
   Aphididae   Aphis   pomi    De Geer 
 Sesiidae    Synanthedon pyri    (Harris)    
 Gracillariidae  Marmara elotella    (Busck) 
 Gracillariidae  Phyllonorycter  crataegella (Clemens)   
  Curculionidae   Anthonomus  quadrigibbus    Say 
 Argyresthiidae  Argyresthia conjugella  Zeller  
 Gracillariidae  Marmara pomonella   Busck   
   Cicadellidae    Empoasca    maligna (Walsh) 
   Pseudococcidae  Phenacoccus aceris  (Signoret)  
 Torymidae   Torymus varians (Walker)    
   Psyllidae   Cacopsylla  mali    (Schmidberger)  
  Curculionidae   Hypothenemus    obscurus    (Fabricius) 
  Bostrichidae    Amphicerus  bicaudatus  (Say)   
 Choreutidae Choreutis   pariana (Clerck)    
 Tischeriidae    Tischeria   malifoliella    Clemens 
   Aphididae   Neophyllaphis   araucariae  Takahashi   
 Argyresthiidae  Argyresthia thuiella    (Packard)   
  Curculionidae   Phyllobius  intrusus    Kono    
 Formicidae  Linepithema humile  (Mayr)  
  Chrysomelidae   Chelymorpha cassidea    (Fabricius) 
  Curculionidae   Ips lecontei    Swaine  
 Noctuidae   Euxoa   auxiliaris  (Grote) 
 Noctuidae   Mythimna    unipuncta   (Haworth)   
 Tenthredinidae  Euura   lasiolepis  Smith   
 Sesiidae    Podosesia   syringae    (Harris)    *Also called lilac borer
  Meloidae    Epicauta    fabricii    (LeConte)   
 Apidae  Bombus  ashtoni (Cresson)   
 Cynipidae   Dryocosmus  kuriphilus  (Yasumatsu) 
   Psyllidae   Diaphorina  citri   Kuwayama    
 Blattellidae    Blattella   asahinai    Mizukubo    
  Cerambycidae    Anoplophora glabripennis    (Motschulsky)   
 Formicidae  Pachycondyla    chinensis   (Emery) 
 Culicidae   Aedes   albopictus  (Skuse) 
  Scarabaeidae    Maladera    castanea    (Arrow) 
 Crambidae   Chilo   suppressalis    (Walker)    
   Diaspididae Aulacaspis  rosarum Borchsenius 
   Aphididae   Brachycorynella asparagi    (Mordvilko) 
  Chrysomelidae   Crioceris   asparagi    (Linnaeus)  
  Gracillariidae  Phyllonorycter  tremuloidiella  (Braun) 
  Chrysomelidae   Chrysomela  crotchi Brown   
   Cicadellidae    Macrosteles quadrcomyza humeralis (
  Pulicidae   Xenopsylla  vexabilis   (Jordan)    
    Ptinidae    Ptinus  ocellus Brown   
 Tetranychidae   Oligonychus punicae (Hirst) 
 Mitranychidae   Oligonychus yothersi    (McGregor)  
   Elachistidae    Stenoma catenifer   Walsingham  
    Thripidae   Scirtothrips    perseae Nakahara    
   Aleyrodidae Trialeurodes    floridensis (Quaintance)eae Comstock    HEMIPTE
 Eotetranychus Pritcharanychus    clitus  Pritchard & Baker   
 Pealius Pealius azaleae (Baker & Moles) 
Thyridopteryx Thyridopteryx   ephemeraeformis (Haworth)   
   Papilionidae    Papilio andraemon bonhotei  Sharpe  
   Pyralidae   Dioryctria  pygmaeella  
   Vespidae    Dolichovespula  maculata    (Linnaeus)  
 Diprionidae Neodiprion  abietis (Hay    COLEOPTE
 Cecidomyiidae   Paradiplosis    tumifex 
 Schizoychidae   Schizotetranychus   celarius    (Banks) 
  Pentalididae    Pentalonia  nigronervosa    Coquerel    
 Cosmonidae  Cosmopolites    sordidus    (Germar)    
 Hesperiidae Erionota    thrax   (Linnaeus)  
 Sesiidae    Podosesia   aureocincta Purrington & Niabrotica balteata LeConte    COLEOPTE
bean fly    Ophiomyia ryon) 
 Noctuidae   Autoplusia  egena   
 Crambidae   Maruca  vitrata (Fabricius) 
 Anthomyiidae    Delia   florilega   (Zetterstedt)   
  Curculionidae   Sternechus  paludatus   (Casey) 
    Thripidae   Caliothrips fasciatus   (Pergande)  
  Chrysomelidae   Acanthoscelides obtectus    (Say)   
   Cicadellidae    Balclutha   saltuella   (Kirschbaum)    
   Cimicidae   Cimex   lectularius Linnaeus    
   Aphididae   Grylloprociphilus   imbricator  (Fitch) 
   Eriococcidae    Cryptococcus    fagisuga    Lindinger   
 Noctuidae   Spodoptera  exigua  
  Chrysomelidae   Erynephala  puncticollis    (Say)   
   Cicadellidae    Circulifer  tenellus    (Baker) 
 Anthomyiidae    Pegomya betae   Curtis  
 Crambidae   Loxostege   sticticalis (Linnaeus)  
 Erebidae    Utetheisa   ornatrix    (Linnaeus)  
 Noctuidae   Mamestra    configurata Walker  
 Tortricidae Epiblema    otiosana    (Clemens)   
 Formicidae  Pheidole    megacephala (Fabricius) 
  Curculionidae   Dryocoetes  betulae Hopkins 
 Coleophoridae   Coleophora  serratella  (Linnaeus)  *Also called cigar casebearer
 Tenthredinidae  Fenusa  pumila  (Leach) 
 Argidae Arge    pectoralis  (Leach) 
 Bucculatricidae Bucculatrix canadensisella  Chambers    
 Pyralidae   Acrobasis   betulella   Hulst   
 Calliphoridae   Protocalliphora     Hough   
 Ixodidae    Haemaphysalis   chordeilis  (Packard)   
 Sphecidae   Sceliphron  caementarium    (Drury) 
 Noctuidae   Actebia fennica (Tauscher)  
  Meloidae    Epicauta    pensylvanica    (De Geer)   
 Calliphoridae   Phormia regina  (Meigen)    
   Thyreocoridae   Corimelaena pulicaria   (Germar)    
 Formicidae  Camponotus  pennsylvanicus  (De Geer)   
  Dermestidae Attagenus   unicolor    (Brahm) 
   Aphididae   Myzus   cerasi  (Fabricius) 
   Aphididae   Toxoptera   aurantii    (Boyer de Fonscolombe)  
 Sphecidae   Dolichurus  stantoni    (Ashmead)   
 Noctuidae   Agrotis ipsilon (Hufnagel)  
  Scarabaeidae    Copris  incertus    Say 
  Curculionidae   Magdalis    barbita (Say)   
    Phlaeothripidae Haplothrips gowdeyi (Franklin)  
  Tenebrionidae   Alphitobius laevigatus  (Fabricius) 
 Cephidae    Trachelus   tabidus (Fabricius) 
 Tabanidae   Tabanus atratus Fabricius   
    Phlaeothripidae Leptothrips mali    (Fitch) 
 Formicidae  Solenopsis  richteri    Forel   
  Dermestidae Dermestes   ater    De Geer 
   Aphididae   Brachycaudus    persicae    (Passerini) 
   Aphididae   Melanocallis    caryaefoliae    (Davis) 
   Diaspididae Nuculaspis  californica (Coleman)   
 Vespidae    Delta   pyriforme philippinense (Bequaert)  
   Coccidae    Saissetia   oleae   (Olivier)   
 Stratiomyidae   Hermetia    illucens    (Linnaeus)  
   Plataspidae Coptosoma   xanthogrammtidae    Apion   occidentale Fall    
 Papilionidae    Papilio polyxenes   Stoll   *Immature called parsleyworm
 Apidae  Bombus  melanopygus Nylander    
   Diaspididae Ischnaspis  longirostris    (Signoret)  
  Scarabaeidae    Ataenius    spretulus   (Haldeman)  
  Curculionidae   Dendroctonus    terebrans   (Olivier)   
  Curculionidae   Xylosandrus compactus   (Eichhoff)  
  Curculionidae   Otiorhynchus    sulcatus    (Fabricius) 
  Curculionidae   Conotrachelus   retentus    (Say)   
 Erebidae    Ascalapha   odorata (Linnaeus)  
  Cleridae    Enoclerus   lecontei    (Wolcott)   
 Heliodinidae    Schreckensteinia    festaliella 
 Lycaenidae  Vaga    blackburni  (Tuely) 
   Nabidae Nabis   blackburni  White   
 Libellulidae    Nesogonia   blackburni  (McLachlan) 
   Cicadellidae    Graminella  nigrifrons  (Forbes)    
 Tenthredinidae  Tethida barda   (Say)   
  Cerambycidae    Callidium   antennatum hesperum lackhorned tree cricket Oecanthus nigricornis F. Walker O
 Vespidae    Vespula consobrina  (Saussure)  
   Ixodidae    Ixodes  scapularis  Say 
  Chrysomelidae   Jonthonota  nigripes    (Olivier)   
   Aphididae   Monellia    caryella    (Fitch) 
 Pyralidae   Dioryctria  clarioralis (Walker)    
   Aphididae   Acyrthosiphon   kondoi  Shinji  
 Pyralidae   Melitara    dentata (Grote)"""

class names(object):
    def __init__(self):
        self.genus = []
        self.species = []
        #f = "names.txt"
        #infile = open(f, 'r')
        lines = namelist.split("\n")
        for line in lines:
            l = line.strip().split()
            self.genus.append(l[1])
            self.species.append(l[2])
    def get_name(self):
        g = random.choice(self.genus)
        s = random.choice(self.species)
        return g+" "+s

if __name__ == "__main__":
    n = names()
    for i in range(10):
        print(n.get_name())
