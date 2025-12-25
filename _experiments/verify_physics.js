
// Mock radialR for Hydrogen (1s)
function radialR(n, l, r, Z = 1) {
    // 1s: R = 2 * Z^(3/2) * exp(-Zr)
    // For Z=1: 2 * exp(-r)
    return 2 * Math.pow(Z, 1.5) * Math.exp(-Z * r);
    // Ignore polynomial parts for simple 1s test
}

const Hydrogen = {
    calculateCumulativePotential: function (n, l, Z, atomType, rMax, steps = 500) {
        const dr = rMax / steps;
        const rValues = new Float32Array(steps);
        const eValues = new Float32Array(steps);

        let integral = 0;
        let prevValue = 0;

        for (let i = 0; i < steps; i++) {
            const r = (i + 1) * dr;
            const R = radialR(n, l, r, Z); // Use simplified mock
            const currentValue = -Z * r * R * R;

            const area = 0.5 * (prevValue + currentValue) * dr;
            integral += area;

            rValues[i] = r;
            eValues[i] = integral;
            prevValue = currentValue;
        }

        return { r: rValues, E: eValues };
    },

    transformHistogramToPotential: function (counts, edges, totalSamples, Z) {
        if (totalSamples === 0) return new Float32Array(counts.length);

        const eValues = new Float32Array(counts.length);
        let cumulativeE = 0;

        for (let i = 0; i < counts.length; i++) {
            const count = counts[i];
            if (count > 0) {
                const r = 0.5 * (edges[i] + edges[i + 1]);
                const prob = count / totalSamples;
                const contribution = prob * (-Z / r);
                cumulativeE += contribution;
            }
            eValues[i] = cumulativeE;
        }
        return eValues;
    }
};

// Test Case 1: Theory
try {
    const res = Hydrogen.calculateCumulativePotential(1, 0, 1, 'H', 5, 100);
    console.log("Theory Calculation Success");
    console.log("Last E:", res.E[99]);
    if (isNaN(res.E[99])) console.error("Theory gave NaN");
} catch (e) {
    console.error("Theory Failed", e);
}

// Test Case 2: Experimental
try {
    const counts = [10, 10, 10];
    const edges = [0, 1, 2, 3];
    const total = 30;
    const expE = Hydrogen.transformHistogramToPotential(counts, edges, total, 1);
    console.log("Exp Calculation Success");
    console.log("Exp Values:", expE);
    if (isNaN(expE[2])) console.error("Exp gave NaN");
} catch (e) {
    console.error("Exp Failed", e);
}
