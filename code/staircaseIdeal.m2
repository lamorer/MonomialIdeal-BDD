buildStaircaseIdeals = (n) -> (
    -- Creamos el anillo
    S := QQ[x_1..x_(4*n+7)];
    
    -- Función auxiliar para crear el producto de variables de forma segura
    makeTerm := (startIdx, len, ringS) -> (
        if len == 0 then return 1;
        if len == 1 then return ringS_(startIdx);
        return product toList apply(0..(len-1), k -> ringS_(startIdx + k));
    );

    offset = 2*n + 3;
    -- Generar I1
    gensI1 := apply(1..n, i -> (
        degX := i + 2;
        degY := 2*n + 5 - degX;
        termX := makeTerm(0, degX, S);          
        termY := makeTerm(offset, degY, S);    
        termX * termY
    ));

    term1 = makeTerm(0,1,S)*makeTerm(offset,2*n + 4,S);
    termLast = makeTerm(0,n + 3,S)*makeTerm(offset,1,S);
    gensI1 = (gensI1,term1,termLast);
    
    
    -- Generar I2
    gensI2 := apply(1..n, i -> (
        degY := i + 1;
        degX := 2*n + 5 - degY;
        termX := makeTerm(0, degX, S);
        termY := makeTerm(offset, degY, S);
        termX * termY
    ));

    term2 = makeTerm(0,2,S)*makeTerm(offset, n + 2, S);
    gensI2 = (gensI2,term2);
    
    return (monomialIdeal gensI1, monomialIdeal gensI2);
);

args = scriptCommandLine;

if #args > 1 then(
    n = value(args#1);
) else (
    n = 10;
)

(I1, I2) = buildStaircaseIdeals(n);

TT = {I1,I2};

out = openOut "BDDs/staircaseIdeals.txt";
out << 2 << endl;
out << 4*n+7 << endl;
for I in TT do{
    out << numgens(I) << endl;
    apply(first entries gens I, i-> (out<<i<<endl));
}
close out;

t0 = cpuTime();
I3 = intersect(I1,I2);
t1 = cpuTime();
tI = t1-t0;
print(I3);


t0 = cpuTime();
I4 = I1+I2;
t1 = cpuTime();
tS = t1-t0;
print(I4);


t0 = cpuTime();
I5 = quotient(I1,I2);
t1 = cpuTime();
tC = t1-t0;
       
print(toString(tI));
print(toString(tS));
print(toString(tC));
