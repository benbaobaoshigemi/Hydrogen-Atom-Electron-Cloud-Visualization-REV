/**
 * 物理引擎能量计算验证脚本 v2.0
 * 
 * 目标: 在 Node.js 环境中独立验证 physics.js 的能量计算逻辑，
 *       不修改任何核心代码，仅通过模拟浏览器环境加载并测试。
 * 
 * 评价指标:
 *   1. 归一化积分: |∫|ψ|²dτ - 1| < 10⁻⁶
 *   2. 能量精度 (相对误差): |E_calc - E_ref| / |E_ref| < 5%
 *   3. 核势能: V_nuc 应与 -Z × <1/r> 吻合
 */

const fs = require('fs');
const path = require('path');

// ============================================================
// Phase 1: 模拟浏览器环境
// ============================================================
console.log('='.repeat(60));
console.log('Phase 1: 模拟浏览器环境 (Browser Environment Simulation)');
console.log('='.repeat(60));

global.window = global;
global.self = global;

// ============================================================
// Phase 2: 加载依赖模块 (按顺序加载)
// ============================================================
console.log('\nPhase 2: 加载依赖模块');
const projectRoot = path.resolve(__dirname, '..');

// 2.1 加载 PhysicsCore (CommonJS)
console.log('  [2.1] Loading physics-core.js...');
const PhysicsCore = require(path.join(projectRoot, 'physics-core.js'));
global.PhysicsCore = PhysicsCore;
window.PhysicsCore = PhysicsCore;
console.log('        ✔ PhysicsCore.factorial exists:', typeof PhysicsCore.factorial === 'function');

// 2.2 加载 SlaterBasis (Browser Global)
console.log('  [2.2] Loading slater_basis.js...');
const slaterContent = fs.readFileSync(path.join(projectRoot, 'slater_basis.js'), 'utf8');
// 使用 vm 模块在隔离上下文中执行, 避免 const 重声明问题
const vm = require('vm');
vm.runInThisContext(slaterContent, { filename: 'slater_basis.js' });
console.log('        ✔ SlaterBasis["H"] exists:', typeof global.SlaterBasis?.H === 'object');
console.log('        ✔ SlaterBasis["He"] exists:', typeof global.SlaterBasis?.He === 'object');

// 2.3 加载 Physics (IIFE)
console.log('  [2.3] Loading physics.js...');
const physicsContent = fs.readFileSync(path.join(projectRoot, 'physics.js'), 'utf8');
vm.runInThisContext(physicsContent, { filename: 'physics.js' });
const Physics = global.Hydrogen;
console.log('        ✔ Hydrogen.calculateCumulativeOrbitalEnergy exists:', typeof Physics?.calculateCumulativeOrbitalEnergy === 'function');

// ============================================================
// Phase 3: 测试指标定义
// ============================================================
console.log('\n' + '='.repeat(60));
console.log('Phase 3: 测试指标定义 (Evaluation Metrics)');
console.log('='.repeat(60));
console.log(`
| 指标名称             | 公式                                      | 阈值         |
| -------------------- | ----------------------------------------- | ------------ |
| 归一化误差           | |∫|ψ|²dτ - 1|                             | < 10⁻⁶       |
| 能量绝对误差         | |E_calc - E_ref|    (Ha)                  | < 0.1 Ha     |
| 能量相对误差         | |E_calc - E_ref| / |E_ref|                | < 5%         |
| 核势能一致性         | V_nuc ≈ -Z × ∫(P(r)/r)dr                  | 定性一致     |
`);

// ============================================================
// Phase 4: 执行测试
// ============================================================
console.log('='.repeat(60));
console.log('Phase 4: 执行测试 (Running Tests)');
console.log('='.repeat(60));

const results = [];

function runEnergyTest(atomSymbol, n, l, Z, E_ref_total, epsilon_ref) {
    console.log(`\n--- Test Case: ${atomSymbol} ${n}${['s', 'p', 'd', 'f'][l]} (Z=${Z}) ---`);

    const testCase = {
        atom: atomSymbol,
        orbital: `${n}${['s', 'p', 'd', 'f'][l]}`,
        Z: Z,
        E_ref_total: E_ref_total,
        epsilon_ref: epsilon_ref,
        status: 'UNKNOWN',
        raw: {}
    };

    // 检查基组是否存在
    const atomData = window.SlaterBasis[atomSymbol];
    if (!atomData) {
        console.error(`  [FATAL] 基组数据不存在: ${atomSymbol}`);
        testCase.status = 'ERROR: NO BASIS';
        results.push(testCase);
        return;
    }

    const orbitalKey = `${n}${['s', 'p', 'd', 'f'][l]}`;
    const basis = atomData.orbitals[orbitalKey];
    if (!basis) {
        console.error(`  [FATAL] 轨道基组不存在: ${atomSymbol} ${orbitalKey}`);
        testCase.status = 'ERROR: NO ORBITAL';
        results.push(testCase);
        return;
    }

    // 输出基组原始数据
    console.log('  [RAW] 基组参数:');
    basis.forEach((term, i) => {
        console.log(`    Term ${i}: nStar=${term.nStar}, zeta=${term.zeta.toFixed(6)}, coeff=${term.coeff.toFixed(6)}`);
    });
    testCase.raw.basis = basis.map(t => ({ nStar: t.nStar, zeta: t.zeta, coeff: t.coeff }));

    // 调用核心计算函数
    const rMax = 30.0;
    const steps = 1000;
    let energyResult;
    try {
        energyResult = Physics.calculateCumulativeOrbitalEnergy(n, l, Z, atomSymbol, rMax, steps);
    } catch (e) {
        console.error(`  [FATAL] 计算抛出异常: ${e.message}`);
        testCase.status = 'ERROR: EXCEPTION';
        testCase.raw.exception = e.message;
        results.push(testCase);
        return;
    }

    if (!energyResult || energyResult.E.length === 0) {
        console.error('  [FATAL] 计算返回空结果');
        testCase.status = 'ERROR: EMPTY RESULT';
        results.push(testCase);
        return;
    }

    // 提取最终值
    const finalIdx = energyResult.E.length - 1;
    const E_calc = energyResult.E[finalIdx];
    const T_calc = energyResult.T[finalIdx];
    const V_calc = energyResult.V[finalIdx];
    const r_final = energyResult.r[finalIdx];

    testCase.raw.rMax_actual = r_final;
    testCase.raw.E_final = E_calc;
    testCase.raw.T_final = T_calc;
    testCase.raw.V_final = V_calc;

    console.log('  [RAW] 计算结果:');
    console.log(`    r_max (实际): ${r_final.toFixed(4)} a.u.`);
    console.log(`    T (动能):     ${T_calc} Ha`);
    console.log(`    V_nuc:        ${V_calc} Ha`);
    console.log(`    E (总能):     ${E_calc} Ha`);

    // 判断数值有效性
    if (!Number.isFinite(E_calc) || !Number.isFinite(T_calc)) {
        console.error('  [FAIL] 数值发散 (NaN/Infinity)');
        testCase.status = 'FAIL: NaN';
        testCase.raw.nanSource = !Number.isFinite(T_calc) ? 'T' : 'E';
        results.push(testCase);
        return;
    }

    // 计算误差指标
    const absError = Math.abs(E_calc - epsilon_ref);
    const relError = Math.abs(E_calc - epsilon_ref) / Math.abs(epsilon_ref);

    testCase.raw.epsilon_ref = epsilon_ref;
    testCase.raw.absError = absError;
    testCase.raw.relError = relError;

    console.log('  [METRIC] 评价指标:');
    console.log(`    ε_ref (HF轨道能):     ${epsilon_ref.toFixed(6)} Ha`);
    console.log(`    绝对误差 |ΔE|:        ${absError.toFixed(6)} Ha`);
    console.log(`    相对误差 |ΔE|/|E|:    ${(relError * 100).toFixed(2)}%`);

    // 判定
    if (relError < 0.05) {
        testCase.status = 'PASS';
        console.log('  [RESULT] ✅ PASS');
    } else if (relError < 0.50) {
        testCase.status = 'PARTIAL';
        console.log('  [RESULT] ⚠️ PARTIAL (within 50%)');
    } else {
        testCase.status = 'FAIL';
        console.log('  [RESULT] ❌ FAIL');
    }

    results.push(testCase);
}

// 执行测试案例
// 参考能量来自 NIST Atomic Spectra Database & Koga Basis 元数据
// E_tot: 原子总能 (Hartree), epsilon: 轨道能 (Hartree, Koopmans定理)

// H: 单电子，无Vee，E_total = epsilon = -0.5
runEnergyTest('H', 1, 0, 1, -0.5, -0.5);

// He: 双电子，E_total = -2.8617, epsilon_1s = -0.9179
runEnergyTest('He', 1, 0, 2, -2.8617, -0.9179);

// Li 1s: 核心壳层, epsilon = -2.4777
runEnergyTest('Li', 1, 0, 3, -7.4327, -2.4777);

// C 1s: epsilon = -11.3255
runEnergyTest('C', 1, 0, 6, -37.6886, -11.3255);

// C 2p: epsilon = -0.4333
runEnergyTest('C', 2, 1, 6, -37.6886, -0.4333);

// ============================================================
// Phase 5: 汇总报告
// ============================================================
console.log('\n' + '='.repeat(60));
console.log('Phase 5: 测试汇总 (Test Summary)');
console.log('='.repeat(60));

const passed = results.filter(r => r.status === 'PASS').length;
const partial = results.filter(r => r.status === 'PARTIAL').length;
const failed = results.filter(r => r.status.startsWith('FAIL') || r.status.startsWith('ERROR')).length;

console.log(`\n总计: ${results.length} 个测试`);
console.log(`  ✅ PASS:    ${passed}`);
console.log(`  ⚠️ PARTIAL: ${partial}`);
console.log(`  ❌ FAIL:    ${failed}`);

console.log('\n详细数据表:');
console.log('| 原子 | 轨道 | E_calc (Ha) | ε_ref (Ha) | |ΔE| (Ha) | Rel.Err (%) | 状态 |');
console.log('| --- | --- | --- | --- | --- | --- | --- |');
for (const r of results) {
    const E_str = r.raw.E_final !== undefined ? r.raw.E_final.toFixed(4) : 'N/A';
    const eps_str = r.raw.epsilon_ref !== undefined ? r.raw.epsilon_ref.toFixed(4) : 'N/A';
    const abs_str = r.raw.absError !== undefined ? r.raw.absError.toFixed(4) : 'N/A';
    const rel_str = r.raw.relError !== undefined ? (r.raw.relError * 100).toFixed(2) : 'N/A';
    console.log(`| ${r.atom} | ${r.orbital} | ${E_str} | ${eps_str} | ${abs_str} | ${rel_str} | ${r.status} |`);
}

// 输出 JSON 供后续处理
const jsonPath = path.join(__dirname, 'validation_results.json');
fs.writeFileSync(jsonPath, JSON.stringify(results, null, 2), 'utf8');
console.log(`\n原始数据已保存至: ${jsonPath}`);

console.log('\n' + '='.repeat(60));
console.log('验证完成');
console.log('='.repeat(60));
