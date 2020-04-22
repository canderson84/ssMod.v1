#simulation prime with A/PR/8/34 and boost with Cal09, test AG Cal09
import random
import myfunc_12Ag_poly4 as func


# get command-line options
vaccine = "monovalent"
#vaccine = "polyvalent"

#immunization parameters
#vacc_interval = 16.0  #days16
#stop_interval = 714.0  #days714
vacc_interval = 180.0  #days180.0
stop_interval = 365.0  #days365.0
ag_initial = 360
#ag_initial = 5
ag_decay = 12.0

#antigen paramaters
stimulation = 1.0
clearance = 0.0001

ep1_clearance = 1.0  #this is relative clearance that modifies the base clearance (set to 0.0001)
ep2_clearance = 1.0
ep3_clearance = 1.0
ep4_clearance = 1.0
ep5_clearance = 1.0
ep6_clearance = 1.0
ep1_immunogenicity = 0.80
ep2_immunogenicity = 1.2
ep3_immunogenicity = 1.2
ep4_immunogenicity = 1.2
ep5_immunogenicity = 1.2
ep6_immunogenicity = 1.2

#immune system parameters
Bcell_initial = 50000000  #see Smith & Perelson PNAS 1999
Bcell_carry = 5000  #see Kuppers et al. EMBO J 1993
Bcell_aff = 10.0
AB_aff = 2.5
tau = 8.0  # time interval for constants (hrs) see Zhang et al Immune Lett. 1988; Liu et al Eur J. Immun. 1991

#####
#####Define Antigens
#AG1 #Cal09
ep1_strain1 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain1 = func.Epitope("poly1", '13431122211221111111', ep2_immunogenicity, ep2_clearance)
ep3_strain1 = func.Epitope("poly2", '24344241232412212342', ep3_immunogenicity, ep3_clearance)
ep4_strain1 = func.Epitope("poly3", '33131232434231143212', ep4_immunogenicity, ep4_clearance)
ep5_strain1 = func.Epitope("poly4", '33333333334323132344', ep5_immunogenicity, ep5_clearance)
ep6_strain1 = func.Epitope("poly5", '14121311132141411114', ep6_immunogenicity, ep6_clearance)
strain1 = func.AntigenType("TypeI", 0)
strain1.add_epitope(ep1_strain1)
strain1.add_epitope(ep2_strain1)
strain1.add_epitope(ep3_strain1)
strain1.add_epitope(ep4_strain1)
strain1.add_epitope(ep5_strain1)
strain1.add_epitope(ep6_strain1)
#AG2 #BR07
ep1_strain2 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain2 = func.Epitope("poly1", '11332114222231111111', ep2_immunogenicity, ep2_clearance)
ep3_strain2 = func.Epitope("poly2", '22122112322332234422', ep2_immunogenicity, ep2_clearance)
ep4_strain2 = func.Epitope("poly3", '31231231241221231131', ep2_immunogenicity, ep2_clearance)
ep5_strain2 = func.Epitope("poly4", '33333313333333333323', ep2_immunogenicity, ep2_clearance)
ep6_strain2 = func.Epitope("poly5", '12121131212121112141', ep2_immunogenicity, ep2_clearance)
strain2 = func.AntigenType("TypeII", 1)
strain2.add_epitope(ep1_strain2)
strain2.add_epitope(ep2_strain2)
strain2.add_epitope(ep3_strain2)
strain2.add_epitope(ep4_strain2)
strain2.add_epitope(ep5_strain2)
strain2.add_epitope(ep6_strain2)
#AG3 #SC18
ep1_strain3 = func.Epitope("con",'44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain3 = func.Epitope("poly1", '13231112211221311111', ep2_immunogenicity, ep2_clearance)
ep3_strain3 = func.Epitope("poly2", '24344241212412241342', ep2_immunogenicity, ep2_clearance)
ep4_strain3 = func.Epitope("poly3", '34131232434131141412', ep2_immunogenicity, ep2_clearance)
ep5_strain3 = func.Epitope("poly4", '33333333331333231333', ep2_immunogenicity, ep2_clearance)
ep6_strain3 = func.Epitope("poly5", '11121311113211112111', ep2_immunogenicity, ep2_clearance)
strain3 = func.AntigenType("TypeIII", 2)
strain3.add_epitope(ep1_strain3)
strain3.add_epitope(ep2_strain3)
strain3.add_epitope(ep3_strain3)
strain3.add_epitope(ep4_strain3)
strain3.add_epitope(ep5_strain3)
strain3.add_epitope(ep6_strain3)
#AG4 #SI06
ep1_strain4 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain4 = func.Epitope("poly1", '11332114223221111111', ep2_immunogenicity, ep2_clearance)
ep3_strain4 = func.Epitope("poly2", '21121112222332233324', ep3_immunogenicity, ep3_clearance)
ep4_strain4 = func.Epitope("poly3", '34231231141221232131', ep4_immunogenicity, ep4_clearance)
ep5_strain4 = func.Epitope("poly4", '31333313333333333323', ep5_immunogenicity, ep5_clearance)
ep6_strain4 = func.Epitope("poly5", '12121141212121112141', ep6_immunogenicity, ep6_clearance)
strain4 = func.AntigenType("TypeIV", 3)
strain4.add_epitope(ep1_strain4)
strain4.add_epitope(ep2_strain4)
strain4.add_epitope(ep3_strain4)
strain4.add_epitope(ep4_strain4)
strain4.add_epitope(ep5_strain4)
strain4.add_epitope(ep6_strain4)
#AG5 #PR34
ep1_strain5 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain5 = func.Epitope("poly1", '11231131224223114211', ep2_immunogenicity, ep2_clearance)
ep3_strain5 = func.Epitope("poly2", '24134144212422214212', ep3_immunogenicity, ep3_clearance)
ep4_strain5 = func.Epitope("poly3", '34332331333131411211', ep4_immunogenicity, ep4_clearance)
ep5_strain5 = func.Epitope("poly4", '32431333334343421333', ep5_immunogenicity, ep5_clearance)
ep6_strain5 = func.Epitope("poly5", '12121121312121122114', ep6_immunogenicity, ep6_clearance)
strain5 = func.AntigenType("TypeV", 4)
strain5.add_epitope(ep1_strain5)
strain5.add_epitope(ep2_strain5)
strain5.add_epitope(ep3_strain5)
strain5.add_epitope(ep4_strain5)
strain5.add_epitope(ep5_strain5)
strain5.add_epitope(ep6_strain5)
#AG6 #NJ76
ep1_strain6 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain6 = func.Epitope("poly1", '13431131211221111111', ep2_immunogenicity, ep2_clearance)
ep3_strain6 = func.Epitope("poly2", '24341241212412214342', ep3_immunogenicity, ep3_clearance)
ep4_strain6 = func.Epitope("poly3", '34132232434131111411', ep4_immunogenicity, ep4_clearance)
ep5_strain6 = func.Epitope("poly4", '33333333331333231333', ep5_immunogenicity, ep5_clearance)
ep6_strain6 = func.Epitope("poly5", '11141311143444444211', ep6_immunogenicity, ep6_clearance)
strain6 = func.AntigenType("TypeVI", 5)
strain6.add_epitope(ep1_strain6)
strain6.add_epitope(ep2_strain6)
strain6.add_epitope(ep3_strain6)
strain6.add_epitope(ep4_strain6)
strain6.add_epitope(ep5_strain6)
strain6.add_epitope(ep6_strain6)
#AG7 #US77
ep1_strain7 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain7 = func.Epitope("poly1", '11341111322221111111', ep2_immunogenicity, ep2_clearance)
ep3_strain7 = func.Epitope("poly2", '21122412212142332221', ep3_immunogenicity, ep3_clearance)
ep4_strain7 = func.Epitope("poly3", '31211213231431231231', ep4_immunogenicity, ep4_clearance)
ep5_strain7 = func.Epitope("poly4", '33233333333333333333', ep5_immunogenicity, ep5_clearance)
ep6_strain7 = func.Epitope("poly5", '12141124242111122111', ep6_immunogenicity, ep6_clearance)
strain7 = func.AntigenType("TypeVII", 6)
strain7.add_epitope(ep1_strain7)
strain7.add_epitope(ep2_strain7)
strain7.add_epitope(ep3_strain7)
strain7.add_epitope(ep4_strain7)
strain7.add_epitope(ep5_strain7)
strain7.add_epitope(ep6_strain7)
#AG8 #BR78
ep1_strain8 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain8 = func.Epitope("poly1", '11341111322221111111', ep2_immunogenicity, ep2_clearance)
ep3_strain8 = func.Epitope("poly2", '21122412212142132221', ep3_immunogenicity, ep3_clearance)
ep4_strain8 = func.Epitope("poly3", '31211213231431231231', ep4_immunogenicity, ep4_clearance)
ep5_strain8 = func.Epitope("poly4", '32333333333333333333', ep5_immunogenicity, ep5_clearance)
ep6_strain8 = func.Epitope("poly5", '12141124242111122111', ep6_immunogenicity, ep6_clearance)
strain8 = func.AntigenType("TypeVIII", 7)
strain8.add_epitope(ep1_strain8)
strain8.add_epitope(ep2_strain8)
strain8.add_epitope(ep3_strain8)
strain8.add_epitope(ep4_strain8)
strain8.add_epitope(ep5_strain8)
strain8.add_epitope(ep6_strain8)
#AG9 #CH83
ep1_strain9 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain9 = func.Epitope("poly1", '11341111322221111111', ep2_immunogenicity, ep2_clearance)
ep3_strain9 = func.Epitope("poly2", '21122412212142132221', ep3_immunogenicity, ep3_clearance)
ep4_strain9 = func.Epitope("poly3", '31211213231431231231', ep4_immunogenicity, ep4_clearance)
ep5_strain9 = func.Epitope("poly4", '33333313333333333333', ep5_immunogenicity, ep5_clearance)
ep6_strain9 = func.Epitope("poly5", '12111122212111122111', ep6_immunogenicity, ep6_clearance)
strain9 = func.AntigenType("TypeIX", 8)
strain9.add_epitope(ep1_strain9)
strain9.add_epitope(ep2_strain9)
strain9.add_epitope(ep3_strain9)
strain9.add_epitope(ep4_strain9)
strain9.add_epitope(ep5_strain9)
strain9.add_epitope(ep6_strain9)
#AG10 #SI86
ep1_strain10 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain10 = func.Epitope("poly1", '11331112222221111111', ep2_immunogenicity, ep2_clearance)
ep3_strain10 = func.Epitope("poly2", '21122112222332233424', ep3_immunogenicity, ep3_clearance)
ep4_strain10 = func.Epitope("poly3", '31211213231431231231', ep4_immunogenicity, ep4_clearance)
ep5_strain10 = func.Epitope("poly4", '33333313333333333333', ep5_immunogenicity, ep5_clearance)
ep6_strain10 = func.Epitope("poly5", '14121141212111122113', ep6_immunogenicity, ep6_clearance)
strain10 = func.AntigenType("TypeX", 9)
strain10.add_epitope(ep1_strain10)
strain10.add_epitope(ep2_strain10)
strain10.add_epitope(ep3_strain10)
strain10.add_epitope(ep4_strain10)
strain10.add_epitope(ep5_strain10)
strain10.add_epitope(ep6_strain10)
#AG11 #BE95
ep1_strain11 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain11 = func.Epitope("poly1", '11323111322121111111', ep2_immunogenicity, ep2_clearance)
ep3_strain11 = func.Epitope("poly2", '21122112222332233424', ep3_immunogenicity, ep3_clearance)
ep4_strain11 = func.Epitope("poly3", '31231231241221231131', ep4_immunogenicity, ep4_clearance)
ep5_strain11 = func.Epitope("poly4", '33333343333333333333', ep5_immunogenicity, ep5_clearance)
ep6_strain11 = func.Epitope("poly5", '12121124242124122111', ep6_immunogenicity, ep6_clearance)
strain11 = func.AntigenType("TypeXI", 10)
strain11.add_epitope(ep1_strain11)
strain11.add_epitope(ep2_strain11)
strain11.add_epitope(ep3_strain11)
strain11.add_epitope(ep4_strain11)
strain11.add_epitope(ep5_strain11)
strain11.add_epitope(ep6_strain11)
#AG12 #NC99
ep1_strain12 = func.Epitope("con", '44444444444444444444', ep1_immunogenicity, ep1_clearance)
ep2_strain12 = func.Epitope("poly1", '11332114222221111114', ep2_immunogenicity, ep2_clearance)
ep3_strain12 = func.Epitope("poly2", '22422112222332233414', ep3_immunogenicity, ep3_clearance)
ep4_strain12 = func.Epitope("poly3", '31231231241221231131', ep4_immunogenicity, ep4_clearance)
ep5_strain12 = func.Epitope("poly4", '33333333333333333333', ep5_immunogenicity, ep5_clearance)
ep6_strain12 = func.Epitope("poly5", '12121121212321122111', ep6_immunogenicity, ep6_clearance)
strain12 = func.AntigenType("TypeXII", 11)
strain12.add_epitope(ep1_strain12)
strain12.add_epitope(ep2_strain12)
strain12.add_epitope(ep3_strain12)
strain12.add_epitope(ep4_strain12)
strain12.add_epitope(ep5_strain12)
strain12.add_epitope(ep6_strain12)

#ADD antigens to list
antigen_list = []
antigen_list.append(strain1)
antigen_list.append(strain2)
antigen_list.append(strain3)
antigen_list.append(strain4)


test_list = []
test_list.append(strain1)
test_list.append(strain2)
test_list.append(strain3)
test_list.append(strain4)
test_list.append(strain5)
test_list.append(strain6)
test_list.append(strain7)
test_list.append(strain8)
test_list.append(strain9)
test_list.append(strain10)
test_list.append(strain11)
test_list.append(strain12)


data_file = 'SC18_CA09_' + str(random.getrandbits(11)) + '.txt'
print "Writing to data file " + data_file

#PRIME SYSTEM
#
#set up antigen population
if (vaccine == "monovalent"):
    ####single serum immunization [initialize]
    #in qoutes changes name in file
    V0 = func.Antigen("V0", 0, strain1)#CA09
    V1 = func.Antigen("V1", 0, strain2)#BR07
    V2 = func.Antigen("V2", ag_initial, strain3)#SC18
    V3 = func.Antigen("V3", 0, strain4)#SI06

if (vaccine == "polyvalent"):
    ###multivalent immunization
    V0 = func.Antigen("V0", 0, strain1)
    V1 = func.Antigen("V1", 0, strain2)
    V2 = func.Antigen("V2", 0, strain3)
    V3 = func.Antigen("V3", 0, strain4)
    V4 = func.Antigen("V4", 0, strain5)
    V5 = func.Antigen("V5", 0, strain6)
    V6 = func.Antigen("V6", 0, strain7)
    V7 = func.Antigen("V7", 0, strain8)
    V8 = func.Antigen("V8", 0, strain9)
    V9 = func.Antigen("V9", 0, strain10)
    V10 = func.Antigen("V10", 0, strain11)
    V11 = func.Antigen("V11", 0, strain12)

#set up B cell populations
nB_gc = func.BCell("GC_B", 0, antigen_list)
nB_stim = func.BCell("Stimulated_B", 0, antigen_list)
nB_me = func.BCell("Memory_B", 0, antigen_list)
nB_pl = func.BCell("Plasma_B", 0, antigen_list)
nAB = func.BCell("Antibody", 0, antigen_list)
nB = func.BCell("Naive_B", Bcell_initial, antigen_list)
nB_llpl = func.BCell("LL_Plasma_B", 0, antigen_list)

GC_population = func.GroupPopulation("GC Bcell Population")
GC_population.add_population(nB_gc)
GC_population.add_population(nB_stim)


#set up population list-all agents in simulation
#stimulated naive cells or memory cells form gc and give rise to plasma cells
PopulationList = []
PopulationList.append(nB_gc)
PopulationList.append(nB_stim)
PopulationList.append(nB_me)
PopulationList.append(nB_pl)
PopulationList.append(nAB)
PopulationList.append(nB)
PopulationList.append(V0)
PopulationList.append(V1)
PopulationList.append(V2)
PopulationList.append(V3)
#PopulationList.append(V4)
#PopulationList.append(V5)
#PopulationList.append(V6)
#PopulationList.append(V7)
#PopulationList.append(V8)
#PopulationList.append(V9)
#PopulationList.append(V10)
#PopulationList.append(V11)

####
####
####
#define equations (rates, counts, time, etc)
#define a system of reactions
A0 = func.TotalReaction()

##Description: Circulating B cell antigen stimulation 
##Naive B cell stimulation is external to the GC
eq1a = func.Stimulation("Free Naive B Cell Stimulation", float(stimulation / 10.0), V0, nB, nB_gc, float(tau / 24.0),
                        Bcell_aff, "immunogenicity")
eq1b = func.Stimulation("Free Naive B Cell Stimulation", float(stimulation / 10.0), V1, nB, nB_gc, float(tau / 24.0),
                        Bcell_aff, "immunogenicity")
eq1c = func.Stimulation("Free Naive B Cell Stimulation", float(stimulation / 10.0), V2, nB, nB_gc, float(tau / 24.0),
                        Bcell_aff, "immunogenicity")
eq1d = func.Stimulation("Free Naive B Cell Stimulation", float(stimulation / 10.0), V3, nB, nB_gc, float(tau / 24.0),
                        Bcell_aff, "immunogenicity")

A0.add_reaction(eq1a)
A0.add_reaction(eq1b)
A0.add_reaction(eq1c)
A0.add_reaction(eq1d)


##Description: Germinal Center B cell antigen stimulation
##             binding affinity with antigen scales exponentially with Hamming Distance
eq2a = func.Stimulation("GC B Cell Stimulation", stimulation, V0, nB_gc, nB_stim, float(tau / 0.25), Bcell_aff,
                        "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2b = func.Stimulation("GC B Cell Stimulation", stimulation, V1, nB_gc, nB_stim, float(tau / 0.25), Bcell_aff,
                        "immunogenicity")
eq2c = func.Stimulation("GC B Cell Stimulation", stimulation, V2, nB_gc, nB_stim, float(tau / 0.25), Bcell_aff,
                        "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2d = func.Stimulation("GC B Cell Stimulation", stimulation, V3, nB_gc, nB_stim, float(tau / 0.25), Bcell_aff,
                        "immunogenicity")
A0.add_reaction(eq2a)
A0.add_reaction(eq2b)
A0.add_reaction(eq2c)
A0.add_reaction(eq2d)

##Description: Naive B cell formation in bone marrow
##NEW 07-20-16
n_antigen = float(len(antigen_list))
multiplier = n_antigen - (n_antigen-1.0)*func.AvgCrossreactivity(antigen_list)
formation_rate = 0.216/multiplier
eq3 = func.Formation("Naive B Cell formation", float(tau / formation_rate), nB)  #produces 500 cells every 108hrs
A0.add_reaction(eq3)
##NEW 07-20-16

##Description: Germinal Center B cell decay, governed by apoptosis. Modeled as a logistic function
##	       with a pre-defined GC population carrying capacity
eq4a = func.PopulationDecay("GC B Cell decay", float(tau / (tau + 1.0)), float(tau / 108.0), Bcell_carry, GC_population,
                            nB_gc)  #base half life of 4.5 days (108hrs)
A0.add_reaction(eq4a)

##Description: Basal B cell decay rate
eq4b = func.Decay("Naive B cell decay", float(tau / 108.0), nB)
A0.add_reaction(eq4b)

##Description: Germinal Center B cell differentiation rate, modeling affinity-dependent T help
eq5 = func.Differentiation("B cell differentiation", 1.0, nB_stim, nB_gc, nB_me, nB_pl, nB_llpl,
                           1.0)  #cell cycle is tau
A0.add_reaction(eq5)

##Description: Antibody production by Plasma cells
eq6a = func.Production("Antibody Production", 1.0, nB_pl, nAB)
eq6b = func.Production("Antibody Production", 1.0, nB_llpl, nAB)
A0.add_reaction(eq6a)
A0.add_reaction(eq6b)

##Description: Antibody decay rate
eq7 = func.Decay("Antibody Decay", float(tau / 360.0), nAB)  #half life of 10 days (360hrs)
A0.add_reaction(eq7)

##Description: Intrinsic and antibody-dependent antigen clearance 
eq8a = func.Decay("Antigen Decay", float(tau / ag_decay), V0)
eq8b = func.Decay("Antigen Decay", float(tau / ag_decay), V1)
eq8c = func.Decay("Antigen Decay", float(tau / ag_decay), V2)
eq8d = func.Decay("Antigen Decay", float(tau / ag_decay), V3)

A0.add_reaction(eq8a)
A0.add_reaction(eq8b)
A0.add_reaction(eq8c)
A0.add_reaction(eq8d)


eq8i = func.Clearance("Antigen Clearance", clearance, V0, nAB, 10000.0, AB_aff, "clearance")
eq8j = func.Clearance("Antigen Clearance", clearance, V1, nAB, 10000.0, AB_aff, "clearance")
eq8k = func.Clearance("Antigen Clearance", clearance, V2, nAB, 10000.0, AB_aff, "clearance")
eq8l = func.Clearance("Antigen Clearance", clearance, V3, nAB, 10000.0, AB_aff, "clearance")

A0.add_reaction(eq8i)
A0.add_reaction(eq8j)
A0.add_reaction(eq8k)
A0.add_reaction(eq8l)


##Description: Plasma cell decay rate
eq9a = func.Decay("Plasma B cell decay", float(tau / 72.0), nB_pl)  #half life of 3 days (72hrs)
eq9b = func.Decay("LL Plasma B cell decay", float(tau / 4800.0), nB_llpl)  #half life of 200 days (4800hrs)
A0.add_reaction(eq9a)
A0.add_reaction(eq9b)

##Description: Circulating memory cell stimulation rate
eq10a = func.Stimulation("Memory B cell Simulation", stimulation, V0, nB_me, nB_stim, float(tau / 24.0), Bcell_aff,
                         "immunogenicity")
eq10b = func.Stimulation("Memory B cell Simulation", stimulation, V1, nB_me, nB_stim, float(tau / 24.0), Bcell_aff,
                         "immunogenicity")
eq10c = func.Stimulation("Memory B cell Simulation", stimulation, V2, nB_me, nB_stim, float(tau / 24.0), Bcell_aff,
                         "immunogenicity")
eq10d = func.Stimulation("Memory B cell Simulation", stimulation, V3, nB_me, nB_stim, float(tau / 24.0), Bcell_aff,
                         "immunogenicity")

A0.add_reaction(eq10a)
A0.add_reaction(eq10b)
A0.add_reaction(eq10c)
A0.add_reaction(eq10d)


#########
#########
#########
#carry out simulation
output = func.FileOutput(data_file, 0.1, PopulationList, antigen_list, test_list)
output.start()

###initial immunization
total_time = float(vacc_interval * 3.0)
t = 0
while ( t <= total_time ):
    dt = A0.MC_TimeStep()
    t += dt
    A0.MC_React()
    output.write(t)

#boost rHA (prime happens at initiation 
if (vaccine == "monovalent"):
    ###single serum
    PopulationList[6].increase(ag_initial)#CA09
if (vaccine == "polyvalent"):
    ##multivalent
    PopulationList[6].increase(ag_initial / 2)#CA09
    PopulationList[7].increase(ag_initial / 2)#BR07

total_time = total_time + float(stop_interval * 3.0)
while ( t <= total_time ):
    dt = A0.MC_TimeStep()
    t += dt
    A0.MC_React()
    output.write(t)

output.finish()
total_time = total_time + float(14.0 * 3.0)
while ( t <= total_time ):
    dt = A0.MC_TimeStep()
    t += dt
    A0.MC_React()

