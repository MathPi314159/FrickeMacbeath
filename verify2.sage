print('This program proves an upper bound for the number of eigenvalues of the Laplacian on the Fricke-Macbeath surface in the interval [0.4, 1.23].\n')

from sage.symbolic.integration.integral import indefinite_integral

a = QQ(0.4)
b = QQ(1.23)
# The endpoints of the interval, stored as rational numbers 
print('a=',a,', b=',b,'\n')

Z = [ QQ(r) for r in [194.67] ]
# The list of double zeros for u, converted to rational numbers
W = [ QQ(r) for r in [1.84, 4.98, 10.02, 16.06, 23.52, 36.65, 47.47] ]
# The list of double zeros for v, converted to rational numbers

print('The list of prescribed double zeros for u is',Z,'\n')
print('The list of prescribed double zeros for v is',W,'\n')

S = (255104361777185173386801383022708600757208600820023523150269595281888727170746101830770724200513723084752076447939620215148426584036900265683951649981818107462106906112699352062859481577229974950638149898343376126809199394203690516618028840393051519562620054808974158837371385861354742353056786219139171089701314927659875897180959884713/197678418206670191567869206336363195993589561606883023205176674291214065322591750888905704197812762927077047423960294759208268950073264960795749316547992113881898555096804191967212957085580555559833271096087716813839240063041326642825872101068814375076288688127852252393681661637431539457817424073181281096241089356063685616661896179124, -85362826800531072617831435894980193282125545089978231597496956539469190542205921019524890163762001171293569948507656795163142333332743569446552520514517837310951475143766820637402750072360302537379764048297948867854052872072284189381590941981968265578736774576263003554488456110045335152397959537608414992061567142506707753359627902769/98839209103335095783934603168181597996794780803441511602588337145607032661295875444452852098906381463538523711980147379604134475036632480397874658273996056940949277548402095983606478542790277779916635548043858406919620031520663321412936050534407187538144344063926126196840830818715769728908712036590640548120544678031842808330948089562, 28499124704363724759023771780308517657947105234076236937626733406800140551850148545337441433798506871469245074039398976658615159146163824028044504546850951692361204362874646079013957759027901816260554148865292978162883145020289420258156751343175306171653379112558426339967080867683597324453184997548830341575733841937400341729303559745/49419604551667547891967301584090798998397390401720755801294168572803516330647937722226426049453190731769261855990073689802067237518316240198937329136998028470474638774201047991803239271395138889958317774021929203459810015760331660706468025267203593769072172031963063098420415409357884864454356018295320274060272339015921404165474044781, -18272124868634385575298731219697534619201385952849962606372895371707935446654634769026674057222744110467877359425621386897946160592795595884930125686716562944096699828225078710833685831530080315023962285754340699203531190063936837085828264025098936815177093142567447798464686190363530337728844598960942576750882685443886294808273096450/49419604551667547891967301584090798998397390401720755801294168572803516330647937722226426049453190731769261855990073689802067237518316240198937329136998028470474638774201047991803239271395138889958317774021929203459810015760331660706468025267203593769072172031963063098420415409357884864454356018295320274060272339015921404165474044781, 78038270846587360573878589683985445667350810837839105228673848076007100533936318987588738165231961230392417882785604563907058959983797181365465815380253543494235397800160185196263861822664341009684311724294848567629300547213773595861936987761027668232297518190502380185544215059705362017460139750368324442500667767187509759365426585000/345937231861672835243771111088635592988781732812045290609059180009624614314535564055584982346172335122384832991930515828614470662628213681392561303958986199293322471419407335942622674899765972229708224418153504424218670110322321624945276176870425156383505204223741441688942907865505194051180492128067241918421906373111449829158318313467, -131697238804481166578496624236951351300513526346054985778111280287271095442574197179029727243783037427812147714358841548552041545758187408956449482695029492598556047542616360353964080879452612152124792366238266508838877938340536384726012022457812259297365604939079167514030891908503398501566829838569624599441457536368446078687983750000/1037811695585018505731313333265906778966345198436135871827177540028873842943606692166754947038517005367154498975791547485843411987884641044177683911876958597879967414258222007827868024699297916689124673254460513272656010330966964874835828530611275469150515612671224325066828723596515582153541476384201725755265719119334349487474954940401, 64491023947123456085105877112837469051473537842491080940600751153518787519559465350398940124410536568572084823540514456739964181025739605589527583104073672756827599235535637018138097089369741540432282517049355098997716933577229667465152449936522250700623678104900422568518824635028725537458815640092150309893082719468567257870812500000/1037811695585018505731313333265906778966345198436135871827177540028873842943606692166754947038517005367154498975791547485843411987884641044177683911876958597879967414258222007827868024699297916689124673254460513272656010330966964874835828530611275469150515612671224325066828723596515582153541476384201725755265719119334349487474954940401, -3227323741526300563718042611286717616160756239026613167887132033604767020272945277059212924786788465287456512865048642523443870345810624027602564749140339706717934834339727337364037191901999244127293574647136907097681004078832033070547813440273727449159968661046229679652147399913347545137390355646105291664207826589634876159375000000/148258813655002643675901904752272396995192171205162267403882505718410548991943813166679278148359572195307785567970221069406201712554948720596811987410994085411423916322603143975409717814185416669874953322065787610379430047280994982119404075801610781307216516095889189295261246228073654593363068054885960822180817017047764212496422134343, -10767191499487257753749334912194034747422109188929890824789994018902910582947280900250825449918907066293613187138560007641957496492784418227848206887688409863691381710173735873421974032405233019381652362258767609215007766564505933487796700667507179769209633312350224654058285574000410173784941020696055016402889660552310000000000000/21179830522143234810843129250324628142170310172166038200554643674058649855991973309525611164051367456472540795424317295629457387507849817228115998201570583630203416617514734853629959687740773809982136188866541087197061435325856426017057725114515825901030930870841312756465892318296236370480438293555137260311545288149680601785203162049, 27579274401054588634065904806260530399891962131190514012826968020426560780811328279014057202599772095793235619781401687715795186009642876589468026126720064592151691994329072805104557755124407790561457148924440761626701957274539455979100235670980911313984456712756368543561984968800046091122705350397969676661724326654500000000000000/2601031818508818310103542188636357842020915284301092410594429924884395596349891459064548739444904775356276939788951246831687749343069275799944069954578843603709191514431634104831749435336586257366227251264312063339990000829491140037182527645642294408898535370103319110443179758387257098129176632541858961792645912579785337061340739199, -622353286798111416847531113531861278340857513593033468801060559029507088464990533000201272114519446969911813181221352320326890897384396617851158123071012666829485929539520689468592537361939921005629820285711664893257589134835981864786730282280984147337276144036901558320012961368828681649175982728529925916307847551656250000000000000/49419604551667547891967301584090798998397390401720755801294168572803516330647937722226426049453190731769261855990073689802067237518316240198937329136998028470474638774201047991803239271395138889958317774021929203459810015760331660706468025267203593769072172031963063098420415409357884864454356018295320274060272339015921404165474044781, 514154359642500267769079736642256108194545644414039425907101621324521017300134718096869761382706919660150763143159554810876695966474114586019493923916485182857712553314269272738027280449711411793295674983970615291477690289801646660915147575865243618091301496193920346600409816226366996954785154571440757675461752239062500000000000000/49419604551667547891967301584090798998397390401720755801294168572803516330647937722226426049453190731769261855990073689802067237518316240198937329136998028470474638774201047991803239271395138889958317774021929203459810015760331660706468025267203593769072172031963063098420415409357884864454356018295320274060272339015921404165474044781, -112477685415091583088822965091498090276499473729101197673609841125702925914130215086605059397308214045059000462745528177057126049389124995354793573609383454835810891789270758793805105311673902903963130175063887570005325132436958215828274912809038753990545949715386240476002170108445175143212066669650468502053740000000000000000000000/16473201517222515963989100528030266332799130133906918600431389524267838776882645907408808683151063577256420618663357896600689079172772080066312443045666009490158212924733682663934413090465046296652772591340643067819936671920110553568822675089067864589690724010654354366140138469785961621484785339431773424686757446338640468055158014927, 57761624100737084452476726714461167187669489989988346226512001906346244009893638488727925865220844573098496738966453711915860053896715150623994005576003021887237725898895643844126475591942772997281361241229912849400134713018997664781582796091573130188423247669979031292235361402626824133392301049588476612869937500000000000000000000/16473201517222515963989100528030266332799130133906918600431389524267838776882645907408808683151063577256420618663357896600689079172772080066312443045666009490158212924733682663934413090465046296652772591340643067819936671920110553568822675089067864589690724010654354366140138469785961621484785339431773424686757446338640468055158014927, -2838445664178907110900049832293765031525135653898958069610183013199616459462312871392400201342646281268251988502918667196263168366561589047328205325419653891448610340661932876164157335178202319017667155668875576087509304086241726250516266230859378744493518946828997289971268758116237299801587339703028731575000000000000000000000000/2353314502460359423427014361147180904685590019129559800061627074895405539554663701058401240450151939608060088380479699514384154167538868580901777577952287070022601846390526094847773298637863756664681798762949009688562381702872936224117525012723980655670103430093479195162876924255137374497826477061681917812393920905520066865022573561, 151809853911372583197583775889682314327332935907494986144962624282275676322425983549956036539267372763622172956090321414805643900621583667965015962982629362997506104274019858760650702877354598762174521578488193456600421028972039314930604874855732107122324980981374811655606939057903813984006331321614418750000000000000000000000000/784438167486786474475671453715726968228530006376519933353875691631801846518221233686133746816717313202686696126826566504794718055846289526967259192650762356674200615463508698282591099545954585554893932920983003229520793900957645408039175004241326885223367810031159731720958974751712458165942159020560639270797973635173355621674191187, 64329518601179143172436587737332735184610243075010594964044105097603491324286160172828458053126006992175876157161998281411520265019771059970084640040633902143914713610934498719740234139007200955894563834954439296472748194272472730267145789617030320391752347919031772717945908694833929448973910267720000000000000000000000000000000/784438167486786474475671453715726968228530006376519933353875691631801846518221233686133746816717313202686696126826566504794718055846289526967259192650762356674200615463508698282591099545954585554893932920983003229520793900957645408039175004241326885223367810031159731720958974751712458165942159020560639270797973635173355621674191187, 583026673888774614953480405490332606278456977814472428152694927861152445593149369722681713244821208958287819696753498740607853271373382194368674977294999143264265252774329792686272856142174101845185286756921537171388828658316668018084824403316951780374406587960311845261830802958578013877880127000000000000000000000000000000000/112062595355255210639381636245103852604075715196645704764839384518828835216888747669447678116673901886098099446689509500684959722263755646709608456092966050953457230780501242611798728506564940793556276131569000461360113414422520772577025000605903840746195401433022818817279853535958922595134594145794377038685424805024765088810598741)
# The vector of coefficients for u (generated with another program) 
print("The vector of coefficients is",S,"\n")

P.<x> = PolynomialRing(QQ) 
# P is the polynomial ring with rational coefficients
L = [ P(gen_laguerre(n,-1/2,x)) for n in range(len(S))]
# A list of the generalized Laguerre polynomials L_n^(-1/2) for n between 0 and len(S)-1

u = sum([S[j]*L[j] for j in range(len(S))])
v = sum([(-1)^j*S[j]*L[j] for j in range(len(S))])
# This is the unique polynomial such that v(x^2)*e^(-x^2/2) is the Fourier transform of u(x^2)*e^(-x^2/2)
print("u is equal to",u,"\n")
print("v is equal to",v,"\n")

u_der = u.derivative()
v_der = v.derivative()
# The first derivative of u and v

u_der_der = u_der.derivative()
v_der_der = v_der.derivative()
# The second derivative of u and v

#The following verifies that u and v vanish at the prescribed orders at the zeros

flag = True
# A Boolean that we modify whenever a condition is not met

for j in range(len(Z)):
    if u(Z[j])!=0 :
        flag = False
        print('u does not vanish at one of the prescribed double zeros.\n')
                        
    if u_der(Z[j])!=0 :
        flag = False
        print('u_der does not vanish at one of the prescribed double zeros.\n')

    if u_der_der(Z[j])>=0 :
        flag = False
        print('either u vanishes to higher order at one of the prescribed zeros or it has a local minimum there, making it positive nearby.\n')

for j in range(len(W)):
    if v(W[j])!=0 :
        flag = False
        print('v does not vanish at one of the prescribed double zeros.\n')
                        
    if v_der(W[j])!=0 :
        flag = False
        print('v_der does not vanish at one of the prescribed double zeros.\n')

    if v_der_der(W[j])<=0 :
        flag = False
        print('either v vanishes to higher order at one of the prescribed zeros or it has a local maximum there, making it negative nearby.\n')

# We then check that u(x^2) has no other zeros in [5.7962,infinity) besides the double ones and is negative at 5.7962.

count = 0
for point in Z:
    if point >= (57962/10000)^2:
        count = count+1


if pari.polsturm(u,[(57962/10000)^2,Infinity])>count: 
    flag = False
    print('u has an extra zero.\n')
    
if u((57962/10000)^2)>=0: 
    flag = False
    print('u(x^2) is non-negative at 5.7962.\n')
    
print('The value of u at 5.7962^2 is',u((57962/10000)^2))

# We also check that that v has no extra zeros in [a-1/4,infinity) and satisfies v(a-1/4)>0.

count = 0
for point in W:
    if point >= a-1/4:
        count = count+1

if pari.polsturm(v,[a-1/4,Infinity])>count:
    flag = False
    print('v has an extra zero.\n')
    
if v(a-1/4)<=0: 
    flag = False
    print('v is non-positive at a-1/4.\n')
                               
# We then check that f^(x) = v(x^2)*e^(-x^2/2) is either monotone or has a unique local maximum on the interval [r(a),r(b)] where r(t)=sqrt(t-1/4)

crit = pari.polsturm(2*v_der-v,[a-1/4,b-1/4])

if crit>1:    
    flag = False
    print('f^ has more than one critical point in [r(a),r(b)].\n')

if crit==1:    
# This means that f^ has a critical point, so we check that it is a local maximum
    
    if 2*v_der(a-1/4)-v(a-1/4)<=0:    
        flag = False
        print('f^ does not have a positive derivative at r(a).\n')
            
    if 2*v_der(b-1/4)-v(b-1/4)>=0:    
        flag = False
        print('f^ does not have a negative derivative at r(b).\n')  

#If the number of critical points is zero, then f^ is monotone on the interval, hence attains its minimum at an endpoint.
        
# After all the checks, we verify the value of flag     
if flag==True:
    print('All the conditions on u and v are satisfied.\n')
    
# Finally, we estimate the upper bound that f yields from above using interval arithmetic
# Here RBF and CBF denote the RealBallField and ComplexBallField from the Arb package

# h is defined so that f^(x) = h(x^2)
def h(t):
    return v(t)*e^(-t/2)

# We first compute the minimum of f^ at r(a) and r(b), which gives a lower bound on the interval [r(a),r(b)].
# For this, we tell Sage to consider these inputs as real ball fields, so that subsequent calculations are done with interval arithmetic. 
lower_bound = min(h(RBF(a-1/4)),h(RBF(b-1/4)))
print('f^ is at least',lower_bound.lower(),'on [r(a),r(b)].\n')

# We then estimate the integral term in the Selberg trace formula
truncated_integral = 12*CBF.integral(lambda x, _: h(x^2)*x*tanh(CBF(pi)*x),0,200)

# For the remainder of the integral, we compute an indefinite integral and verify that it is a polynomial multiplied by e^(-x^2/2)
V(x) = indefinite_integral(h(x^2)*x,x)
D = e^(x^2/2)*V(x)
p = D.full_simplify()
if p.is_polynomial(x)!=True:
    print('Error, V is not of the form e^(-x^2/2)p(x) where p is a polynomial.\n')
remainder = -12*tanh(pi*200)*V(RBF(200))

# We compute the contribution of the shortest and second closed geodesics in the geometric term 
def f(x):
    return u(x^2)*e^(-x^2/2)
def w(x): 
    return (x/sinh(x/2))*f(x)
# This is the contribution of each geodesic of length x in the geometric sum (the factor of 2 is not there because we consider unoriented geodesics)

LenC=arccosh(cot(pi/3)*cot(pi/7))
LenB=arcsinh(sin(pi/3)*sinh(LenC))
LenA=arcsinh(sin(pi/7)*sinh(LenC))
l1 = 4*(LenA+LenB+LenC) 
# This is the systole of the Fricke-Macbeath surface 
l2 = 4*arccosh(sinh(LenA+LenB+LenC)*sinh(LenA+LenB+LenC))
# This is the length of the second shortest geodesics on the Fricke-Macbeath surface 
geom = RBF((1/sqrt(2*pi))*(126*w(l1)+252*w(l2)))
# There are 126 systoles of length l1 and 252 curves of length l2, so geom is the resulting sum

# We then compute the resulting upper bound on the number of eigenvalues in [a,b] (recall that h(-1/4) = f^(i/2))
multiplicity = (truncated_integral + remainder + geom - h(RBF(-1/4)))/lower_bound
print('The final upper bound is at most',RBF(multiplicity).upper())
# We convert the CBF element multiplicity to RBF to get an upper bound on it 
# (only CBF can compute integrals, but its balls do not have upper bounds)
