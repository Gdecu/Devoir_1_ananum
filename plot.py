import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def import_data(filename):
    # Charger le fichier en ignorant les lignes d'en-tête répétées
    data = []
    with open(filename, "r") as f:
        next(f)  # Ignorer la première ligne d'en-tête
        for line in f:
            if not line.startswith("n"):  # Ignorer les éventuelles répétitions d'en-tête
                data.append(list(map(float, line.strip().split(","))))

    # Convertir en tableau NumPy
    array = np.array(data)
    print(array)    
    return data


"""### Plot : k constant et n varie
data = import_data("times_k.txt")
data = np.array(data)
print(data)
n = data[:, 0]
nx = data[:, 1]
ny = data[:, 2] # = k
t = data[:, 3]   # en secondes
print(n, nx, ny, t)

# Plot taille en fonction du temps
# Application du style seaborn
sns.set(style="darkgrid")  # Style quadrillé similaire
plt.figure(figsize=(8, 6))
plt.plot(n, t, marker='o', linestyle='-', color='r')  # en normal
#plt.plot(np.log10(n), np.log10(t), marker='o', linestyle='-', color='b')  # en log-log base 10
#plt.xlabel("Taille de la matrice : log(n)", fontsize=12)
plt.xlabel("Taille de la matrice : log(n)", fontsize=12)
#plt.ylabel("Temps d'exécution log(ms)", fontsize=12)
plt.ylabel("Temps d'exécution s", fontsize=12)
plt.title("Complexité temporelle de qr_eigs avec k = {}".format(ny[0]), fontsize=14)
plt.grid(True, linestyle="--", linewidth=1, alpha=0.7)
plt.tight_layout()
plt.savefig("complexite_temporelle.png")"""


### Plot : n constant et k varie
data = import_data("times_n.txt")
data = np.array(data)
print(data)
k = data[:, 2]  # k varie
t = data[:, 3]
n_constant = data[0, 0]  # n est constant

# Plot k en fonction du temps
sns.set(style="darkgrid")  # Application du style seaborn
plt.figure(figsize=(8, 6))
plt.plot(k, t, marker='o', linestyle='-', color='g')  # en normal
plt.xlabel("Valeur de k", fontsize=12)
plt.ylabel("Temps d'exécution (s)", fontsize=12)
plt.title("Complexité temporelle de qr_eigs avec n = {}".format(n_constant), fontsize=14)
plt.grid(True, linestyle="--", linewidth=1, alpha=0.7)
plt.tight_layout()
plt.savefig("complexite_temporelle_k_varie.png")




"""

# Plot taille en fonction du temps
plt.figure(figsize=(8, 6))
plt.plot(x, y, marker='o', linestyle='-', color='r', label="Données")
plt.xlabel("Taille de la matrice", fontsize=12)
plt.ylabel("Temps d'exécution (s)", fontsize=12)
plt.title("Complexité temporelle de qr_eigs", fontsize=14)
plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("complexite_temporelle.png")

plt.clf()  # Clear the current figure


# Log-log plot avec échelle logarithmique
plt.figure(figsize=(8, 6))
plt.loglog(x, y, marker='o', linestyle='-', color='b', label="Données")
plt.xlabel("Taille de la matrice (log10)", fontsize=12)
plt.ylabel("Temps d'exécution (s) (log10)", fontsize=12)
plt.title("Complexité temporelle de qr_eigs (échelle logarithmique)", fontsize=14)
plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig("log_complexite_temporelle.png")"""