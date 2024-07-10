print('This program proves a lower bound for the number of eigenvalues of the Laplacian on the Fricke-Macbeath surface in the interval [0.77,1.28].\n')

from sage.symbolic.integration.integral import indefinite_integral

a = QQ(0.77)
b = QQ(1.28)
# The endpoints of the interval, stored as rational numbers 
print('a=',a,', b=',b,'\n')

Z = [ QQ(r) for r in [] ]
# The list of double zeros for u, converted to rational numbers
W = [ QQ(r) for r in [1.99, 4.99, 9.92, 16.54, 23.28, 32.14, 46.56, 56.0] ]
# The list of double zeros for v, converted to rational numbers

print('The list of prescribed double zeros for u is',Z,'\n')
print('The list of prescribed double zeros for v is',W,'\n')

S = (9989031108052300674714523734481548094342607332138164649813412019056123961874163983/3358343884094078149981113481664652767531909516157919120618754460166596276258741696, -11526373396012652686515098766944717113185894340746579897335589095360273330683469/4883895124938506035139112678253725600588445885475556627023572059999016266373104, 5156946667326841286175817351448693521902613986504852611047395565404705658877629/2817402587327246770118383793342829502963011339058656980384861124300835802230488, -13369906820551403325062721742013292388749203295603958649762354835845559795676700/9462546804568355443082030814116593351140113903108890964858170866248094016097889, 177804440255941911483159104035915143530971590235117740672496281771819110354989375/164918672879619909150858251331746341262727699454183528244670977954609638566277494, -469497055358799733062599115817963152254797540020154990354111578258696573010156250/577215355078669682028003879661112194419546948089642348856348422841133734981971229, 116489261607332437548666139761541181737042019676495208547278799614327693093750000/192405118359556560676001293220370731473182316029880782952116140947044578327323743, -12186875302110381522653532560633856326270388091335121779351056333898859375000000/27486445479936651525143041888624390210454616575697254707445162992434939761046249, 8733754246603273022631095042377582109699487211480140596327131627006562500000000/27486445479936651525143041888624390210454616575697254707445162992434939761046249, -224922576595880268039828496368390264397545690488030701854956738093750000000000/1018016499256913019449742292171273711498319132433231655831302333053145917075787, 149705757886249139660637816631732096781565078903717979780876812500000000000000/1018016499256913019449742292171273711498319132433231655831302333053145917075787, -2429793210735338740832650726160923603661648256973758485937500000000000000/26374167705300992757577716836479538627900182192109423970343851733287026013, 4923454043489829853100517383828276357277065034525597171601562500000000000000/92546954477901183586340208379206701045301739312111968711936575732104174279617, -2556642388786663855842200606076807832006343636587735550781250000000000000000/92546954477901183586340208379206701045301739312111968711936575732104174279617, 164868972258434176873051595933095133590444855739289062500000000000000000000/13220993496843026226620029768458100149328819901730281244562367961729167754231, -1856744776013172084434840341753358760995814299218750000000000000000000000/400636166571000794746061508135093943919055148537281249835223271567550538007, 5648648401969174124076582213887789709849085625000000000000000000000000000/4406997832281008742206676589486033383109606633910093748187455987243055918077, -54644949230649209539571618111315260588562500000000000000000000000000000/259235166604765220129804505263884316653506272582946691069850352190767995181)
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

#The following verifies that u and v vanish at the precribed orders at the zeros

flag = True
# A Boolean that we modify whenever a condition is not met

for j in range(len(Z)):
    if u(Z[j])!=0 :
        flag = False
        print('u does not vanish at one of the prescribed double zeros.\n')
                        
    if u_der(Z[j])!=0 :
        flag = False
        print('u_der does not vanish at one of the prescribed double zeros.\n')

    if u_der_der(Z[j])<=0 :
        flag = False
        print('either u vanishes to higher order at one of the prescribed zeros or it has a local maximum there, making it negative nearby.\n')

for j in range(len(W)):
    if v(W[j])!=0 :
        flag = False
        print('v does not vanish at one of the prescribed double zeros.\n')
                        
    if v_der(W[j])!=0 :
        flag = False
        print('v_der does not vanish at one of the prescribed double zeros.\n')

    if v_der_der(W[j])>=0 :
        flag = False
        print('either v vanishes to higher order at one of the prescribed zeros or it has a local minimum there, making it positive nearby.\n')

# We then check that u is non-negative on [0,infinity).

if pari.polsturm(u,[0,Infinity])>len(Z): 
    flag = False
    print('u has an extra zero.\n')
    
if u(0)<=0: 
    flag = False
    print('u is non-positive at the origin.\n')

# We also check that that v has no extra zeros in [b-1/4,infinity) and satisfies v(b-1/4)<0. Note that b-1/4 is smaller than any of the prescribed double roots.

if pari.polsturm(v,[b-1/4,Infinity])>len(W):
    flag = False
    print('v has an extra zero.\n')
    
if v(b-1/4)>=0: 
    flag = False
    print('v is non-negative b-1/4.\n')
        
print('The value of v at b-1/4 is',v(b-1/4))
  
                               
# We then check that f^(x) = v(x^2)*e^(-x^2/2) is decreasing on the interval [r(a),r(b)] where r(t)=sqrt(t-1/4)

if pari.polsturm(2*v_der-v,[a-1/4,b-1/4])>0:    
    flag = False
    print('f^ has a critical point in [r(a),r(b)].\n')

if 2*v_der(a-1/4)-v(a-1/4)>=0:    
    flag = False
    print('f^ does not have a negative derivative at r(a).\n')

# After all the checks, we verify the value of flag     
if flag==True:
    print('All the conditions on u and v are satisfied.\n')
    
# Finally, we estimate the lower bound that f yields using interval arithmetic
# Here RBF and CBF denote the RealBallField and ComplexBallField from the Arb package

# h is defined so that f^(x) = h(x^2)
def h(t):
    return v(t)*e^(-t/2)

# Since f^ is decreasing on [r(a),r(b)], the value f^(r(a)) is an upper bound on that interval. We compute its value using interval arithmetic.
upper_bound = h(RBF(a-1/4))
print('f^ is at most',upper_bound.upper(),'on [r(a),r(b)].\n')

# We then estimate the integral term in the Selberg trace formula
truncated_integral = 12*CBF.integral(lambda x, _: h(x^2)*x*tanh(CBF(pi)*x),0,200)

# For the remainder of the integral, we compute an indefinite integral and verify that it is a polynomial multiplied by e^(-x^2/2)
V(x) = indefinite_integral(h(x^2)*x,x)
D = e^(x^2/2)*V(x)
p = D.full_simplify()
if p.is_polynomial(x)!=True:
    print('Error, V is not of the form e^(-x^2/2)p(x) where p is a polynomial.\n')
remainder = -12*V(RBF(200))

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

# We then compute the resulting lower bound on the number of eigenvalues in [a,b] (recall that h(-1/4) = f^(i/2))
multiplicity = (truncated_integral + remainder + geom - h(RBF(-1/4)))/upper_bound
print('The final lower bound is at least',RBF(multiplicity).lower())
# We convert the CBF element multiplicity to RBF to get a lower bound on it 
# (only CBF can compute integrals, but its balls do not have lower bounds)
