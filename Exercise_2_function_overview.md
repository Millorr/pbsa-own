# Übersicht der Funktionen in Prog2_1Simulation, Prog2_2Simulation und Prog2_3Simulation

## Gemeinsame Funktionen in allen Versionen

### Mesh- und Cloth-Building
- **`buildMesh()`**:
  - Erzeugt ein Gitter von Vertices und Indizes für die Visualisierung.

- **`buildCloth()`**:
  - Initialisiert Positionen, Geschwindigkeiten und Kräfte der Partikel im Gitter.

### Rendering und Darstellung
- **`resize(int w, int h)`**:
  - Aktualisiert die Projektionsmatrix basierend auf dem Seitenverhältnis.

- **`render()`**:
  - Führt die Simulation aus und rendert das Gitter.
  - Bindet den Shader und setzt Uniforms (wie `modelView` und `projection`).

### Simulation
- **`step()`**:
  - Führt einen Simulationsschritt durch (inkl. Kraftberechnung, Integration und Randbedingungen).

- **`computeForces()`**:
  - Berechnet Gravitations- und Federkräfte.

- **`applyBoundaryConditions()`**:
  - Fixiert Ecken des Gitters gemäß den Randbedingungen.

### Hilfsfunktionen
- **`getIndex(int i, int j)`**:
  - Wandelt 2D-Indizes in 1D-Indizes um.

- **`getVectorIndex(int i, int j)`**:
  - Wandelt 2D-Indizes in 1D-Vektor-Indizes um.

### Reset
- **`reset(double dt0Param, BoundaryCondition boundaryConditionParam, ...)`**:
  - Setzt die Simulation mit neuen Parametern zurück.

---

## Unterschiede und neue Funktionen

### Unterschiede in `Aufgabe b)`
- **`step()`**:
  - Hinzugefügt: Sparse-Matrix-Konstruktionen und Jacobi-Matrizen für Kräfte.
  - Verwendet einen linearen Solver, um Geschwindigkeitsänderungen (`delta_v`) zu berechnen.

- **`addSpringForceAndJacobians(...)`** *(Neu)*:
  - Berechnet Kräfte und Jacobi-Matrizen für ein Federpaar und fügt diese zum Gesamtsystem hinzu.

- **`computeForcesAndJacobians(...)`** *(Neu)*:
  - Berechnet alle Kräfte und Jacobi-Matrizen für das Gitter.

### Unterschiede in `Aufgabe c)`
- **`render()`**:
  - Simulationsschritte werden mehrfach pro Frame ausgeführt (`maxSimSteps`).
  - Echtzeitfaktor (`realTimeFactor`) wird berechnet.

- **`step()`**:
  - Übernimmt die Struktur von `Prog2_1Simulation`, jedoch mit Anpassungen für mehrfaches Ausführen pro Frame.

- **`reset(...)`**:
  - Zusätzlicher Parameter `maxSimSteps`.

---

## Zusammenfassung der Änderungen
- In `b)`:
  - Erweiterung um Jacobian-basierte Kräfteberechnung.
  - Hinzufügen von Jacobi-Matrizen und linearem Solver.
- In `c)`:
  - Anpassungen für mehrfache Simulationsschritte pro Frame.
  - Hinzufügen von `maxSimSteps`-Logik.
