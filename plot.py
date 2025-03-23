import numpy as np
import matplotlib.pyplot as plt

def import_data(filename="times.txt"):
    data = np.loadtxt(filename)
    return data


# Plot la complexité temporelle en fct du temps d'exécution
data = import_data("times.txt")[0:]
n = data.shape[0]
x = (np.arange(n)+1)**2
y = data

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
plt.savefig("log_complexite_temporelle.png")