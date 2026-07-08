import subprocess
import pandas as pd #para ejecutar, activar entorno virtual: (pyhton3 -m venv venv), source venv/bin/activate (pip install matplotlib)
import matplotlib.pyplot as plt


def plotCSV(t=-1, archivo="BDDs/timesIntersection.csv", ejeX="ratio", titulo=None):
    if t == -1:
        cols = ["ratio", "variables", "gens II", "gens I1", "gens I2", "CoCoALib", "BDD"]
        colores = {
            "CoCoALib": "blue",
            "BDD": "green"
        }
        tiempos = ["CoCoALib", "BDD"]
    else:
        cols = ["ratio", "variables", "gens II", "gens I1", "gens I2", "CoCoALib", "BDD", "M2"]
        colores = {
            "CoCoALib": "blue",
            "BDD": "green",
            "M2": "red"
        }
        tiempos = ["CoCoALib", "BDD", "M2"]
    
    try:
        data = pd.read_csv(archivo, usecols=cols)
    except ValueError as e:
        print(f"Error al leer columnas: {e}")
        print("Intentando leer todas las columnas...")
        data = pd.read_csv(archivo)
        # Verificar que tenemos las columnas necesarias
        required_cols = ["variables"] + tiempos
        for col in required_cols:
            if col not in data.columns:
                print(f"Error: Columna '{col}' no encontrada en el archivo")
                print(f"Columnas disponibles: {list(data.columns)}")
                return


    min_var_total = int(data["variables"].min())
    max_var_total = int(data["variables"].max())
    
    
    rangos = []
    for inicio in range(min_var_total, max_var_total + 1, 20):
        fin = inicio + 19
        rangos.append((inicio,fin))
    
    if rangos and rangos[-1][1] > max_var_total:
        rangos[-1] = (rangos[-1][0], max_var_total)
    
    
    for var_min, var_max in rangos:
        subset = data[(data["variables"] >= var_min) & (data["variables"] <= var_max)]
        
        if len(subset) == 0:
            print(f"No hay datos para variables {var_min}-{var_max}")
            continue
        
        subset = subset.sort_values(ejeX)
        
        line_styles = ['-', '--', ':']
        marker_styles = ['o', 'v', 'X']

        fig, ax = plt.subplots(figsize=(10, 6))

        for i, tiempo in enumerate(tiempos):
            ax.plot(
                subset[ejeX],
                subset[tiempo],
                marker=marker_styles[i % len(marker_styles)],
                color=colores[tiempo],
                linestyle=line_styles[i % len(marker_styles)],
                linewidth=1.5,
                markersize=4,
                label=tiempo
            )

        ax.set_xlabel(ejeX)
        ax.set_ylabel("tiempo (s)")
        if titulo:
            title_text = f"{titulo}\nVariables {var_min}-{var_max} (n={len(subset)})"
        else:
            title_text = f"Variables {var_min}-{var_max}: {ejeX} vs Tiempo (n={len(subset)})"
        
        ax.set_title(title_text, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=14)

        plt.yscale("log")
        plt.xscale("log")
        plt.tight_layout()
        plt.show()


def _configurar_grafico(ax, subset, ejeX, tiempos, colores, titulo, info_extra=""):
    """Función interna para aplicar el formato visual a los gráficos."""
    line_styles = ['-', '--', ':']
    marker_styles = ['o', 'v', 'X']

    for i, tiempo in enumerate(tiempos):
        ax.plot(
            subset[ejeX],
            subset[tiempo],
            marker=marker_styles[i % len(marker_styles)],
            color=colores.get(tiempo, "black"),
            linestyle=line_styles[i % len(line_styles)],
            linewidth=1.5,
            markersize=4,
            label=tiempo
        )

    ax.set_xlabel(fr"$|G(I)|+|G(I^\vee)|$", fontsize=16)
    ax.set_ylabel("time (s)",fontsize = 16)
    ax.set_title(f"{titulo}", fontweight='bold', fontsize=22)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=14)
    ax.set_yscale("log")
    ax.set_xscale("log")


def plotCSV_completo(t=-1, archivo="BDDs/timesIntersection.csv", ejeX="ratio", titulo="Colon ideal times"):
    """Dibuja todos los datos del CSV en un solo gráfico."""
    data, tiempos, colores = _preparar_datos(t, archivo)
    if data is None: return

    data = data.sort_values(ejeX)
    fig, ax = plt.subplots(figsize=(12, 7))
    
    var_min, var_max = data["variables"].min(), data["variables"].max()
    info = f"Rango total: {var_min}-{var_max} variables"
    
    _configurar_grafico(ax, data, ejeX, tiempos, colores, titulo, info)
    plt.tight_layout()
    plt.show()

def plotCSV_por_tramos(t=-1, archivo="BDDs/timesIntersection.csv", ejeX="ratio", titulo="Gráfico", salto=20):
    """Dibuja los datos divididos en ventanas según el número de variables."""
    data, tiempos, colores = _preparar_datos(t, archivo)
    if data is None: return

    min_var = int(data["variables"].min())
    max_var = int(data["variables"].max())
    
    for inicio in range(min_var, max_var + 1, salto):
        fin = inicio + salto - 1
        subset = data[(data["variables"] >= inicio) & (data["variables"] <= fin)]
        
        if len(subset) == 0: continue
        
        subset = subset.sort_values(ejeX)
        fig, ax = plt.subplots(figsize=(10, 6))
        
        info = f"Variables {inicio}-{fin}"
        _configurar_grafico(ax, subset, ejeX, tiempos, colores, titulo, info)
        
        plt.tight_layout()
        plt.show()

def _preparar_datos(t, archivo):
    """Carga el CSV y define columnas según el parámetro t."""
    if t == -1:
        cols = ["ratio", "variables", "gens II", "gens I1", "gens I2", "CoCoALib", "BDD","Sum_gen"]
        tiempos = ["CoCoALib", "BDD"]
    else:
        cols = ["ratio", "variables", "gens I_dual", "gens I1", "CoCoALib", "BDD", "M2","Sum_gen"]
        tiempos = ["CoCoALib", "BDD", "M2"]
    
    colores = {"CoCoALib": "blue", "BDD": "green", "M2": "red"}

    try:
        data = pd.read_csv(archivo, sep=';',usecols=lambda x: x in cols)
        return data, tiempos, colores
    except Exception as e:
        print(f"Error: {e}")
        return None, None, None






# Compilar

resultado1 = subprocess.run("make BDDs/monIdeal_operations", shell=True, capture_output=True, text=True)
print("=== MAKE OUTPUT ===")
print(resultado1.stdout)
print(resultado1.stderr)

for j in range (1,10):
    for i in range (2,15):
        # Ejecutar Macaulay2
        s = "m2 --script BDDs/randomMonIdealSF.m2 y " + str(j) + " " + str(i)
        resultado3 = subprocess.run(s, shell=True, capture_output=True, text=True)
        print("=== M2 OUTPUT ===")
        print(resultado3.stdout)
        print(resultado3.stderr)

        lines = resultado3.stdout.strip().split('\n')
        if len(lines) >= 4:
            tI, tS, tC, tD, tP = float(lines[0]), float(lines[1]), float(lines[2]), float(lines[3]), float(lines[4])
        else:
            tI, tS, tC, tD, tP = -1, -1, -1, -1
        
        s = "BDDs/monIdeal_operations " + str(tI) + " " + str(tS) + " " + str(tC) + " " + str(tD) + " " + str(tP)
        resultado2 = subprocess.run(s, shell=True, capture_output=True, text=True)
        print("=== EXECUTION OUTPUT ===")
        print(resultado2.stdout)
        print(resultado2.stderr)


#for i in range (1,100):
        # Ejecutar Macaulay2
        #s = "m2 --script BDDs/staircaseIdeal.m2 " + str(i)
        #resultado3 = subprocess.run(s, shell=True, capture_output=True, text=True)
        #print("=== M2 OUTPUT ===")
        #print(resultado3.stdout)
        #print(resultado3.stderr)

        #lines = resultado3.stdout.strip().split('\n')
        #if len(lines) >= 3:
        #    tI, tS, tC = float(lines[0]), float(lines[1]), float(lines[2])
        #else:
        #    tI, tS, tC = -1, -1, -1, -1, -1
        
        #s = "BDDs/monIdeal_operations " + str(tI) + " " + str(tS) + " " + str(tC)
        #resultado2 = subprocess.run(s, shell=True, capture_output=True, text=True)
        #print("=== EXECUTION OUTPUT ===")
        #print(resultado2.stdout)
        #print(resultado2.stderr)


#plotCSV_completo(tI,"BDDs/timesIntersection_staircase.csv","ratio","Intersection times")
#plotCSV_completo(tS,"BDDs/timesSum_staircase.csv","ratio","Sum times")
#plotCSV_completo(tC,"BDDs/timesColon_staircase.csv","ratio", "Colon ideal times")

#plotCSV_completo(tI, "BDDs/timesIntersection.csv", "ratio", "Tiempos Intersección (Global)")

#plotCSV_completo(0, "BDDs/timesDual.csv", "Sum_gen", "Alexander Dual times")