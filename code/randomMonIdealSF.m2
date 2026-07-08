--polarización de ideales (no square-free) aleatorios
--m2 --script BDDs/randomMonIdealSF.m2
loadPackage("RandomMonomialIdeals");

args = scriptCommandLine;

if #args > 3 then(
    n = value(args#2);
    d = value(args#3);
) else (
    n = 2;
    d = 3;
)

N=2;

createExponentialProbVector = (d, minProb, maxProb) -> (
    apply(toList(1..d), deg -> (
        t = (deg-1)/(d-1);
        baseProb = minProb * ((maxProb/minProb)^t);
        noiseFactor = 0.1 + 0.3 * t;
        perturbation = random(1.0 - noiseFactor, 1.0 + noiseFactor);
        max(minProb, min(maxProb, baseProb * perturbation))
    ))
);

minProb = 0.01;
maxProb = 0.7;

pVector = createExponentialProbVector(d, minProb, maxProb);

found = false;
while not found do (
    II = randomMonomialIdeals(n,d,pVector,N); 
    found = (II_0 != 0 and II_1 != 0);
)

JJ=apply(II,i->polarize(i));

--JJ es una lista de ideales squarefree, CADA UNO EN UN ANILLO DISTINTO

T=QQ[x_1..x_(n*d)]

--T es un anillo grande (el número de variables hay que sacarlo de II (o de JJ))

TT = {};
for I in JJ do{
    soporte=apply(first entries gens I,j->apply(support(j),i->index i));
    Tsoporte=apply(soporte,j->apply(j,i->T_i));
    TT = append(TT,monomialIdeal(apply (Tsoporte,i->product i)));
}

out = openOut "BDDs/randomIdeal.txt";
out << N << endl;
out << n*d << endl;
for I in TT do{
    out << numgens(I) << endl;
    apply(first entries gens I, i-> (out<<i<<endl));
}

close out;

if #args > 1 then (
    opt = args#1;
    if opt=="y" then(
	I1 = TT_0;
	I2 = TT_1;
	t0 = cpuTime();
	--I3 = intersect(I1,I2);
	t1 = cpuTime();
	tI = t1-t0;

	t0 = cpuTime();
	--I4 = I1+I2;
	t1 = cpuTime();
	tS = t1-t0;

	t0 = cpuTime();
	--I5 = quotient(I1,I2);
	t1 = cpuTime();
	tC = t1-t0;

    t0 = cpuTime();
    I6 = dual I1;
    t1 = cpuTime();
    tD = t1-t0;

    t0 = cpuTime();
    primaryDecomposition I1;
    t1 = cpuTime();
    tP = t1 - t0;
        
	print(toString(tI));
	print(toString(tS));
	print(toString(tC));
    print(toString(tD));
    print(toString(tP));
    ) else (
	print("-1","-1","-1","-1","-1");
    )
) else (
   print("-1","-1","-1","-1","-1");
)
