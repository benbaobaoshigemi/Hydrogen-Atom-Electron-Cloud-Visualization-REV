// Hydrogen atom physics model (real orbitals). Attaches to window.Hydrogen
(function(){
  const A0 = 1; // Bohr radius unit scale for visualization
  const PI = Math.PI;
  const TWO_PI = 2*Math.PI;

  const FACT = (()=>{ const f=[1]; for(let i=1;i<=64;i++) f[i]=f[i-1]*i; return f; })();
  function factorial(n){ return FACT[n] ?? Infinity; }
  function binomialInt(n,k){ if(k<0||k>n) return 0; return factorial(n)/ (factorial(k)*factorial(n-k)); }

  // Generalized Laguerre L_k^{alpha}(x), alpha integer >=0, small k
  function generalizedLaguerre(k, alpha, x){
    let sum=0;
    for(let i=0;i<=k;i++){
      const c = ((i%2)?-1:1) * binomialInt(k+alpha, k-i) / factorial(i);
      sum += c * Math.pow(x,i);
    }
    return sum;
  }

  // Associated Legendre P_l^m(x), Condon–Shortley
  function associatedLegendre(l, m, x){
    const mm = Math.abs(m);
    if(l<mm) return 0;
    let pmm = 1.0;
    if(mm>0){
      const s = Math.sqrt(Math.max(0,1-x*x));
      let fact = 1.0;
      for(let i=1;i<=mm;i++){ pmm *= -fact * s; fact += 2.0; }
    }
    if(l===mm) return pmm;
    let pmmp1 = x * (2*mm+1) * pmm;
    if(l===mm+1) return pmmp1;
    let pll = 0;
    for(let ll=mm+2; ll<=l; ll++){
      pll = ((2*ll-1)*x*pmmp1 - (ll+mm-1)*pmm)/(ll-mm);
      pmm = pmmp1; pmmp1 = pll;
    }
    return pll;
  }

  // Complex spherical harmonic Y_l^m(θ,φ) with Condon–Shortley phase
  function Ylm_complex(l, m, theta, phi){
    const mm = Math.abs(m);
    const Plm = associatedLegendre(l, mm, Math.cos(theta));
    const N = Math.sqrt(((2*l+1)/(4*PI)) * (factorial(l-mm)/factorial(l+mm)));
    const base = N * Plm;
    if(m===0){
      return { re: base, im: 0 };
    }
    const cos_mphi = Math.cos(mm*phi);
    const sin_mphi = Math.sin(mm*phi);
    // e^{imφ} = cos(mφ) + i sin(mφ)
    // apply Condon-Shortley phase: included in associatedLegendre sign already
    // For m<0, Y_l^{-m} = (-1)^m (Y_l^m)^*
    if(m>0){
      return { re: base * cos_mphi, im: base * sin_mphi };
    } else {
      const sign = (mm % 2) ? -1 : 1;
      // (-1)^m with m negative equals (-1)^{|m|}
      return { re: sign * base * cos_mphi, im: -sign * base * sin_mphi };
    }
  }

  function Ylm_abs2(l, m, theta, phi){
    const y = Ylm_complex(l, m, theta, phi);
    return y.re*y.re + y.im*y.im;
  }

  // Real spherical harmonics (linear combos of complex):
  // m=0: Y_l0 real
  // m>0: cos-type: √2 Re Y_l^m,  sin-type: √2 Im Y_l^m
  function realYlm_abs2(l, m, type, theta, phi){
    if(m===0){
      return Ylm_abs2(l, 0, theta, phi);
    }
    const mm = Math.abs(m);
    const y = Ylm_complex(l, mm, theta, phi);
    if(type==='c'){
      const v = Math.SQRT2 * y.re; // √2 Re
      return v*v;
    } else { // 's'
      const v = Math.SQRT2 * y.im; // √2 Im
      return v*v;
    }
  }

  // Real spherical harmonics value (not squared), normalized so ∫|Y|^2 dΩ = 1
  function realYlm_value(l, m, type, theta, phi){
    if(m===0){
      // Y_l^0 is real
      const y = Ylm_complex(l, 0, theta, phi);
      return y.re;
    }
    const mm = Math.abs(m);
    const y = Ylm_complex(l, mm, theta, phi);
    if(type==='c'){
      return Math.SQRT2 * y.re;
    } else {
      return Math.SQRT2 * y.im;
    }
  }

  // Normalized radial R_nl(r) with ∫|R|^2 r^2 dr = 1
  function radialR(n,l,r,Z=1,a0=A0){
    if(n<=0||l<0||l>=n) return 0;
    const rho = (2*Z*r)/(n*a0);
    const k = n - l - 1;
    if(k<0) return 0;
    const pref = Math.pow(2*Z/(n*a0), 1.5) * Math.sqrt( factorial(n-l-1) / (2*n*factorial(n+l)) );
    const poly = generalizedLaguerre(k, 2*l+1, rho);
    return pref * Math.exp(-rho/2) * Math.pow(rho, l) * poly;
  }

  function radialPDF(n,l,r,Z=1,a0=A0){
    const R = radialR(n,l,r,Z,a0);
    return r*r*(R*R);
  }

  function density3D_real(angKey, n,l, r, theta, phi, Z=1,a0=A0){
    const R = radialR(n,l,r,Z,a0);
    let Y2 = 1/(4*PI);
    if (angKey && typeof angKey === 'object'){
      Y2 = realYlm_abs2(angKey.l, angKey.m, angKey.t, theta, phi);
    }
    return (R*R)*Y2;
  }

  function recommendRmax(n, a0=A0){ return 15*n*n*a0; }

  function radialGrid(n,l,rmax,num=512,Z=1,a0=A0){
    const rs = new Float32Array(num);
    const ps = new Float32Array(num);
    const dr = rmax/(num-1);
    for(let i=0;i<num;i++){
      const r = i*dr; rs[i]=r; ps[i]=radialPDF(n,l,r,Z,a0);
    }
    return { r: rs, Pr: ps };
  }

  function histogramRadialFromSamples(rArray, bins=80, rmax=null, normalize=true){
    const N = rArray.length; if(N===0) return {edges:[],counts:[]};
    const maxr = rmax ?? Math.max(...rArray);
    
    // 确保最小范围和bins数
    const effectiveMaxr = Math.max(maxr, 0.1);
    const effectiveBins = Math.max(bins, 10);
    
    const edges = new Float32Array(effectiveBins+1);
    const counts = new Float32Array(effectiveBins);
    const dr = effectiveMaxr/effectiveBins;
    
    for(let i=0;i<=effectiveBins;i++) edges[i]=i*dr;
    for(let i=0;i<N;i++){
      const r = rArray[i]; 
      if(r<0||r>effectiveMaxr) continue;
      const b = Math.min(effectiveBins-1, Math.floor(r/dr)); 
      counts[b]+=1;
    }
    
    if(normalize){
      let area=0; 
      for(let i=0;i<effectiveBins;i++) area += counts[i]*dr;
      if(area>0){ 
        for(let i=0;i<effectiveBins;i++) counts[i]/=area; 
      }
    }
    return { edges, counts, dr, rmax:effectiveMaxr };
  }

  function histogramThetaFromSamples(thetaArray, bins=90, normalize=true){
    const N = thetaArray.length; if(N===0) return {edges:[],counts:[]};
    const edges = new Float32Array(bins+1);
    const counts = new Float32Array(bins);
    const dth = Math.PI/bins;
    for(let i=0;i<=bins;i++) edges[i]=i*dth;
    for(let i=0;i<N;i++){
      const t = thetaArray[i]; if(t<0||t>Math.PI) continue;
      const b = Math.min(bins-1, Math.floor(t/dth)); counts[b]+=1;
    }
    if(normalize){ let area=0; for(let i=0;i<bins;i++) area += counts[i]*dth; if(area>0){ for(let i=0;i<bins;i++) counts[i]/=area; } }
    return { edges, counts, dθ:dth };
  }

  function orbitalParamsFromKey(key){
    // Map UI key to {n,l,angKey:{l,m,t}}
    const R = (n,l,m,t)=>({n,l,angKey:{l,m,t}});
    switch(key){
      // n=1
      case '1s': return R(1,0,0,'c');
      // n=2
      case '2s': return R(2,0,0,'c');
      case '2pz': return R(2,1,0,'c');
      case '2px': return R(2,1,1,'c');
      case '2py': return R(2,1,1,'s');
      // n=3
      case '3s': return R(3,0,0,'c');
      case '3pz': return R(3,1,0,'c');
      case '3px': return R(3,1,1,'c');
      case '3py': return R(3,1,1,'s');
      case '3d_z2': return R(3,2,0,'c');
      case '3d_xz': return R(3,2,1,'c');
      case '3d_yz': return R(3,2,1,'s');
      case '3d_xy': return R(3,2,2,'s');
      case '3d_x2-y2': return R(3,2,2,'c');
      // n=4
      case '4s': return R(4,0,0,'c');
      case '4pz': return R(4,1,0,'c');
      case '4px': return R(4,1,1,'c');
      case '4py': return R(4,1,1,'s');
      case '4d_z2': return R(4,2,0,'c');
      case '4d_xz': return R(4,2,1,'c');
      case '4d_yz': return R(4,2,1,'s');
      case '4d_xy': return R(4,2,2,'s');
      case '4d_x2-y2': return R(4,2,2,'c');
      // seven real f
      case '4f_z3': return R(4,3,0,'c');
      case '4f_xz2': return R(4,3,1,'c');
      case '4f_yz2': return R(4,3,1,'s');
      case '4f_z(x2-y2)': return R(4,3,2,'c');
      case '4f_xyz': return R(4,3,2,'s');
      case '4f_x(x2-3y2)': return R(4,3,3,'c');
      case '4f_y(3x2-y2)': return R(4,3,3,'s');
      default: return R(1,0,0,'c');
    }
  }

  function orbitalKey(params) {
    return `${params.n}${params.l}${params.angKey.m}${params.angKey.t}`;
  }

  // 计算角向边缘概率密度 P(θ) = sin(θ) × ∫|Y|² dφ
  // 对于实球谐函数，无论 m=0 还是 m≠0，公式统一为 2π × N² × P_lm² × sin(θ)
  function angularPDF_Theta(l, m, theta) {
    const mm = Math.abs(m);
    const Plm = associatedLegendre(l, mm, Math.cos(theta));
    const N2 = ((2*l+1)/(4*PI)) * (factorial(l-mm)/factorial(l+mm));
    return 2 * PI * N2 * Plm * Plm * Math.sin(theta);
  }

  window.Hydrogen = {
    A0,
    radialR,
    radialPDF,
    density3D_real,
  realYlm_abs2,
  realYlm_value,
    radialGrid,
    histogramRadialFromSamples,
    histogramThetaFromSamples,
    recommendRmax,
    orbitalParamsFromKey,
    orbitalKey,
    angularPDF_Theta, // Export new function
  };
})();
