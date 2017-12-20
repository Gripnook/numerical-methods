import csv

nodes = dict()
potentials = dict()


# Computes the energy in the square with (x, y) as the bottom-left corner.
def compute_energy(x, y):
    epsilon = 8.85e-12
    v1 = potentials[nodes[(x + 1, y)]]
    v2 = potentials[nodes[(x, y)]]
    v3 = potentials[nodes[(x, y + 1)]]
    v4 = potentials[nodes[(x + 1, y + 1)]]
    V = [v1, v2, v3, v4]
    S = [[1.0, -0.5, 0.0, -0.5],
         [-0.5, 1.0, -0.5, 0.0],
         [0.0, -0.5, 1.0, -0.5],
         [-0.5, 0.0, -0.5, 1.0]]
    # Energy = 0.5 * e0 * V^T S V
    energy = 0.0
    for i in range(len(S)):
        for j in range(len(S[0])):
            energy += V[i] * S[i][j] * V[j]
    return 0.5 * epsilon * energy


if __name__ == '__main__':
    # Fill node dictionary.
    index = 1
    for x in range(6):
        for y in range(6):
            if x >= 4 and y == 5:
                continue
            nodes[(x, y)] = index
            index += 1

    # Get potentials.
    with open('report/question-2/results.dat', 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            potentials[int(row[0])] = float(row[3])

    # Compute energy.
    energy = 0.0
    for x in range(5):
        for y in range(5):
            if x >= 3 and y == 4:
                continue
            energy += compute_energy(x, y)
    energy *= 4

    # Compute capacitance.
    # E = 0.5 * C * V^2 => C = 2 * E / V^2
    capacitance = 2 * energy / 15 ** 2
    print capacitance
