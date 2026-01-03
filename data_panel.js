// 数据面板逻辑（与主渲染逻辑解耦）
// 暴露 DataPanel 全局对象，供 script.js 调用
(function () {
  console.log('[data_panel.js] build 2026-01-03 v11 (EnergyCumulativeDiag nearCore tol)');
  const state = {
    chart: null,
    collapsed: false,
    live: false, // 由主程序控制是否实时更新
    isResized: false, // 标记是否被用户手动调整过大小
    waveHidden: true, // 记录波函数曲线的隐藏状态（默认隐藏）
    angularTheoryHidden: false, // 记录角向理论曲线的隐藏状态
    compareTheoryHidden: true, // 记录比照模式理论曲线的隐藏状态（默认隐藏）
    phiTheoryHidden: false, // 记录φ角向理论曲线的隐藏状态
    chart: null,                // Chart.js 实例
    // 滚动生成相关缓存：用于稳定坐标轴/网格，避免每次滚动导致“硬刷新式跳变”
    rollingCache: {
      lastEnabled: null,
      // key -> {min, max}
      yAxisLocks: {},
      // type -> {rmax, bins}
      compareRGridLocks: {}
    }
  };

  // 统一调整所有曲线线宽：在“原来的 1/3”基础上整体放大 1.5 倍 => 1/2
  const LINE_WIDTH_SCALE = (1 / 3) * 1.5;

  function applyThinLineStyle(datasets) {
    if (!Array.isArray(datasets)) return;
    for (const ds of datasets) {
      if (!ds) continue;
      // Chart.js 在折线图中 dataset.type 通常为空（默认继承 chart.type='line'）。
      // 为了让 compare 模式也能统一细线样式，这里把“看起来像折线”的 dataset 也纳入。
      const data0 = Array.isArray(ds.data) && ds.data.length > 0 ? ds.data[0] : null;
      const looksLikeXY = !!(data0 && typeof data0 === 'object' && ('x' in data0) && ('y' in data0));
      const isLine = (ds.type === 'line') || looksLikeXY || (ds.showLine === true) || (ds.showLine !== false && (ds.borderDash || ds.tension !== undefined));
      if (!isLine) continue;
      if (typeof ds.borderWidth === 'number') {
        ds.borderWidth = ds.borderWidth * LINE_WIDTH_SCALE;
      }
    }
  }

  function isRollingEnabled() {
    return !!(window.ElectronCloud && window.ElectronCloud.state && window.ElectronCloud.state.rollingMode && window.ElectronCloud.state.rollingMode.enabled);
  }

  function syncRollingCache() {
    const enabled = isRollingEnabled();
    const cache = state.rollingCache;
    if (cache.lastEnabled === null) {
      cache.lastEnabled = enabled;
      return enabled;
    }
    if (cache.lastEnabled && !enabled) {
      cache.yAxisLocks = {};
      cache.compareRGridLocks = {};
    }
    cache.lastEnabled = enabled;
    return enabled;
  }

  function computeYRangeFromDatasets(datasets) {
    let min = Infinity;
    let max = -Infinity;
    if (!Array.isArray(datasets)) return null;

    for (const ds of datasets) {
      if (!ds || !Array.isArray(ds.data)) continue;
      for (const point of ds.data) {
        let y = null;
        if (typeof point === 'number') {
          y = point;
        } else if (point && typeof point === 'object' && typeof point.y === 'number') {
          y = point.y;
        }
        if (!Number.isFinite(y)) continue;
        if (y < min) min = y;
        if (y > max) max = y;
      }
    }

    if (!Number.isFinite(min) || !Number.isFinite(max)) return null;
    if (min === max) {
      // 避免 min=max 导致 Chart.js 退化
      const pad = (min === 0) ? 1 : Math.abs(min) * 0.05;
      return { min: min - pad, max: max + pad };
    }
    const pad = (max - min) * 0.05;
    return { min: min - pad, max: max + pad };
  }

  function applyRollingAnimationPolicy(chart) {
    if (!chart || !chart.options) return;
    const rolling = isRollingEnabled();
    if (!chart.options.animation) chart.options.animation = {};
    // 滚动模式：完全关闭动画，避免肉眼“硬刷新”
    chart.options.animation.duration = rolling ? 0 : (chart.options.animation.duration ?? 300);
  }

  function applyRollingYAxisLock(chart, datasets, lockKey, typeHint) {
    if (!chart || !chart.options || !chart.options.scales || !chart.options.scales.y) return;

    const rolling = isRollingEnabled();
    if (!rolling) return;

    // 角向/φ 已有固定 0-1，不再二次锁定
    if (typeHint === 'angular' || typeHint === 'phi' || typeHint === 'azimuthal') return;

    const cache = state.rollingCache;
    const existing = cache.yAxisLocks[lockKey];
    const yRange = existing || computeYRangeFromDatasets(datasets);
    if (!yRange) return;

    // 概率类图：强制从 0 起，防止 x 轴（基线）被“抬升”
    // 仅锁 y.max，y.min 固定为 0
    const isProbabilityLike = (typeHint === 'radial');
    if (isProbabilityLike) {
      const locked = existing || { min: 0, max: yRange.max };
      if (!existing) cache.yAxisLocks[lockKey] = locked;
      chart.options.scales.y.min = 0;
      chart.options.scales.y.max = locked.max;
      return;
    }

    if (!existing) cache.yAxisLocks[lockKey] = yRange;
    chart.options.scales.y.min = yRange.min;
    chart.options.scales.y.max = yRange.max;
  }

  function init() {
    // === 数据面板滑出逻辑 ===
    const dataPanel = document.getElementById('data-panel');
    const dataCollapseBtn = document.getElementById('data-collapse-btn');
    const dataTab = document.getElementById('data-panel-tab');
    const enlargeBtn = document.getElementById('chart-enlarge-btn');
    const resizeHandle = document.getElementById('chart-resize-handle');
    const chartContainer = dataPanel ? dataPanel.querySelector('.chart-container') : null;

    // 固定间距常量（与CSS变量--spacing-panel一致）
    const PANEL_SPACING = 16;
    const MIN_CHART_HEIGHT = 220;

    // 计算最小面板高度
    function calculateMinPanelHeight() {
      if (!dataPanel || !chartContainer) return 400;
      const chartTopOffset = chartContainer.getBoundingClientRect().top - dataPanel.getBoundingClientRect().top;
      // 最小面板高度 = 图表顶部偏移 + 最小图表高度 + 底部padding(16px)
      return chartTopOffset + MIN_CHART_HEIGHT + PANEL_SPACING;
    }

    // 延迟计算以确保CSS已加载
    let minPanelHeight = 400;
    requestAnimationFrame(() => {
      minPanelHeight = Math.max(calculateMinPanelHeight(), dataPanel ? dataPanel.getBoundingClientRect().height : 400);
    });

    if (dataCollapseBtn && dataPanel) {
      dataCollapseBtn.addEventListener('click', () => {
        state.collapsed = true;
        dataPanel.classList.add('collapsed');
      });
    }

    if (dataTab && dataPanel) {
      dataTab.addEventListener('click', () => {
        state.collapsed = false;
        dataPanel.classList.remove('collapsed');
      });
    }

    // 复位面板大小的函数
    function resetPanelSize() {
      // 清除所有手动设置的样式
      dataPanel.style.width = '';
      dataPanel.style.height = '';

      // 清除图表容器的手动样式
      if (chartContainer) {
        chartContainer.style.height = '';
      }

      state.isResized = false;

      // 恢复图标
      if (enlargeBtn) {
        enlargeBtn.textContent = '⤢';
        enlargeBtn.title = '放大/还原图表';
      }
    }

    function updateChartContainerHeight(panelHeight) {
      if (!chartContainer || !dataPanel) {
        return;
      }
      const panelContent = dataPanel.querySelector('.panel-content');
      const panelBody = dataPanel.querySelector('.panel-body');
      if (!panelContent || !panelBody) return;

      // 使用固定的16px间距（与CSS变量--spacing-panel一致）
      const PANEL_SPACING = 16;

      // 获取面板边框
      const panelStyle = window.getComputedStyle(dataPanel);
      const panelBorderTop = parseFloat(panelStyle.borderTopWidth) || 0;
      const panelBorderBottom = parseFloat(panelStyle.borderBottomWidth) || 0;

      // 计算面板头部占用的高度
      const panelHeader = dataPanel.querySelector('.panel-header');
      let usedHeight = 0;

      if (panelHeader) {
        usedHeight += panelHeader.getBoundingClientRect().height;
        usedHeight += PANEL_SPACING; // header的margin-bottom
      }

      // 计算其他控件占用的高度（不包含图表容器）
      // 数据面板中有2个control-group，每个后面有16px间距
      const controlGroups = panelBody.querySelectorAll('.control-group');
      const controlGroupCount = controlGroups.length;
      controlGroups.forEach((group, index) => {
        usedHeight += group.getBoundingClientRect().height;
        // 每个control-group后面都有margin-bottom（16px）
        // 因为图表容器不是control-group，所以最后一个control-group也有margin
        usedHeight += PANEL_SPACING;
      });

      // 图表容器的margin-top已经包含在上面最后一个control-group的margin-bottom中
      // 所以不需要再加

      // 计算图表容器可用高度
      // 面板内容区域 = panelHeight - 边框 - 上下padding(各16px)
      // 图表高度 = 内容区域 - usedHeight
      const contentAreaHeight = panelHeight - panelBorderTop - panelBorderBottom - PANEL_SPACING * 2;
      const availableHeight = contentAreaHeight - usedHeight;
      const safeHeight = Math.max(MIN_CHART_HEIGHT, availableHeight);
      chartContainer.style.height = `${safeHeight}px`;
    }

    if (enlargeBtn && dataPanel) {
      enlargeBtn.addEventListener('click', () => {
        // 如果处于手动调整大小状态，点击按钮则复位
        if (state.isResized) {
          resetPanelSize();
          // 确保仍在 enlarged 模式
          if (!dataPanel.classList.contains('enlarged')) {
            dataPanel.classList.add('enlarged');
          }
        } else {
          // 正常的切换逻辑
          const wasEnlarged = dataPanel.classList.contains('enlarged');
          dataPanel.classList.toggle('enlarged');

          // 退出 enlarged 模式时也要复位
          if (wasEnlarged) {
            resetPanelSize();
          }
        }

        // 更新按钮图标：放大时显示向内箭头，缩小时显示向外箭头
        const isEnlarged = dataPanel.classList.contains('enlarged');
        enlargeBtn.textContent = isEnlarged ? '⤡' : '⤢';
        enlargeBtn.title = isEnlarged ? '还原图表' : '放大图表';
      });
    }

    // === 拖动调整大小逻辑 ===
    if (resizeHandle && dataPanel && enlargeBtn && chartContainer) {
      let isResizing = false;
      let offsetX = 0;
      let offsetY = 0;

      resizeHandle.addEventListener('mousedown', (e) => {
        if (!dataPanel.classList.contains('enlarged')) return;

        isResizing = true;
        const handleRect = resizeHandle.getBoundingClientRect();
        offsetX = e.clientX - handleRect.left;
        offsetY = e.clientY - handleRect.top;

        dataPanel.classList.add('resizing');
        e.preventDefault();
      });

      document.addEventListener('mousemove', (e) => {
        if (!isResizing) return;

        const handleTargetX = e.clientX - offsetX;
        const handleTargetY = e.clientY - offsetY;
        // 抓手在面板外24px处，所以面板顶部位置比抓手位置低24px
        const panelTopY = handleTargetY + 24;
        const newWidth = Math.max(300, window.innerWidth - handleTargetX - 20 - 24);
        const rawHeight = window.innerHeight - panelTopY - 20;
        const newHeight = Math.max(minPanelHeight, rawHeight);

        dataPanel.style.width = `${newWidth}px`;
        dataPanel.style.height = `${newHeight}px`;
        updateChartContainerHeight(newHeight);

        if (!state.isResized) {
          state.isResized = true;
          enlargeBtn.textContent = '⟲';
          enlargeBtn.title = '复位大小';
        }
      });

      document.addEventListener('mouseup', () => {
        if (!isResizing) return;
        isResizing = false;
        dataPanel.classList.remove('resizing');
        window.dispatchEvent(new Event('resize'));
      });
    }

    // === 控制面板滑出逻辑 ===
    const controlPanel = document.getElementById('control-panel');
    const controlCollapseBtn = document.getElementById('control-collapse-btn');
    const controlTab = document.getElementById('control-panel-tab');

    if (controlCollapseBtn && controlPanel) {
      controlCollapseBtn.addEventListener('click', () => {
        controlPanel.classList.add('collapsed');
      });
    }

    if (controlTab && controlPanel) {
      controlTab.addEventListener('click', () => {
        controlPanel.classList.remove('collapsed');
      });
    }
  }

  function ensureChart() {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) {
      return null;
    }
    if (state.chart) {
      return state.chart;
    }

    // 强制设置canvas背景为透明
    ctx.style.backgroundColor = 'transparent';

    const rolling = isRollingEnabled();
    state.chart = new Chart(ctx, {
      type: 'bar',
      data: { labels: [], datasets: [] },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        backgroundColor: 'transparent',
        animation: {
          duration: rolling ? 0 : 300,
          easing: 'easeOutQuart'
        },
        interaction: {
          intersect: false,
          mode: 'index'
        },
        layout: {
          padding: 0
        },
        plugins: {
          legend: {
            labels: {
              // color: '#e8e8e8', // 移除硬编码颜色，改由 generateLabels 动态控制
              font: { size: 11, weight: '500' },
              padding: 12,
              generateLabels: function (chart) {
                const original = Chart.defaults.plugins.legend.labels.generateLabels;
                const labels = original.call(this, chart);
                labels.forEach(label => {
                  if (label.hidden) {
                    label.fontColor = '#666666'; // 非激活时变灰
                    label.fillStyle = '#666666';
                    label.strokeStyle = '#666666';
                    label.textDecoration = 'line-through';
                  } else {
                    label.fontColor = '#e8e8e8'; // 激活时亮色
                  }
                });
                return labels;
              }
            }
          },
          tooltip: {
            backgroundColor: 'rgba(20, 20, 20, 0.95)',
            titleColor: '#ffffff',
            bodyColor: '#e8e8e8',
            borderColor: 'rgba(255, 255, 255, 0.2)',
            borderWidth: 1,
            cornerRadius: 6
          }
        },
        scales: {
          x: {
            ticks: {
              color: '#d0d0d0',
              font: { size: 10 },
              maxTicksLimit: 12 // 限制标签数量，避免拥挤
            },
            grid: {
              color: 'rgba(255,255,255,0.08)',
              lineWidth: 1
            },
            border: {
              color: 'rgba(255,255,255,0.15)'
            }
          },
          y: {
            ticks: {
              color: '#d0d0d0',
              font: { size: 10 }
            },
            grid: {
              color: 'rgba(255,255,255,0.08)',
              lineWidth: 1
            },
            border: {
              color: 'rgba(255,255,255,0.15)'
            }
          }
        }
      }
    });
    console.log('新图表实例创建成功'); // 调试信息
    return state.chart;
  }

  function renderChartRadial(hist, theory) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) {
      return;
    }

    // 在更新或重建图表前，先从现有图表实例中捕获可见性状态
    if (state.chart && state.chart.data && state.chart.data.datasets) {
      const waveDatasetIndex = state.chart.data.datasets.findIndex(ds => ds.label === '波函数（幅值）');
      if (waveDatasetIndex !== -1) {
        state.waveHidden = !state.chart.isDatasetVisible(waveDatasetIndex);
      }

    }

    // 若当前是极坐标图，需销毁后重新创建柱状图
    if (state.chart && state.chart.config && state.chart.config.type !== 'bar') {
      try { state.chart.destroy(); } catch (e) { console.warn('图表销毁失败:', e); }
      state.chart = null;
    }
    const chart = ensureChart();
    if (!chart) {
      return;
    }
    // 将直方图中心作为 x 轴（线性 r），理论曲线在相同中心取样
    const centers = (theory && theory.centers) || (() => { const n = hist.counts.length; const a = new Array(n); for (let i = 0; i < n; i++) a[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]); return a; })();

    // 动态调整显示精度
    const maxValue = Math.max(...centers);
    const decimalPlaces = maxValue > 100 ? 1 : (maxValue > 10 ? 2 : 3);

    chart.data.labels = centers.map(v => v.toFixed(decimalPlaces));
    const datasets = [
      {
        label: '径向概率密度 (归一化)',
        data: Array.from(hist.counts),
        backgroundColor: 'rgba(255,255,255,0.7)',
        borderColor: 'rgba(255,255,255,0.95)',
        borderWidth: 0,
        barPercentage: 1.0,
        categoryPercentage: 1.0,
        borderRadius: 2,
        borderSkipped: false,
        order: 10,
      },
    ];
    if (theory && theory.values && theory.values.length) {
      datasets.push({
        type: 'line',
        label: '理论曲线',
        data: theory.values.map((y, i) => ({ x: i, y })),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y',
        tension: 0.2,
        order: 5,
      });
    }
    applyThinLineStyle(datasets);
    chart.data.datasets = datasets;

    // 添加坐标轴配置
    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = {
          display: true,
          text: '径向距离 (a₀)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = {
          display: true,
          text: '概率密度',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };

        // 【关键修复】从角向/φ图切换过来时，清除固定的 0-1 轴范围
        delete chart.options.scales.y.min;
        delete chart.options.scales.y.max;
      }

      // 【关键修复】清除能量图残留的 y2 轴
      delete chart.options.scales.y2;
    }

    // 优化x轴显示
    if (chart.options.scales && chart.options.scales.x && chart.options.scales.x.ticks) {
      chart.options.scales.x.ticks.autoSkip = true;
      chart.options.scales.x.ticks.maxTicksLimit = 15;
    }

    applyRollingAnimationPolicy(chart);
    applyRollingYAxisLock(chart, datasets, 'bar:radial', 'radial');

    try {
      chart.update('none');
    } catch (error) {
      console.error('图表更新失败:', error);
      state.chart.destroy();
      state.chart = null;
      const newChart = ensureChart();
      if (newChart) {
        newChart.data.datasets = datasets;
        applyRollingAnimationPolicy(newChart);
        applyRollingYAxisLock(newChart, datasets, 'bar:radial', 'radial');
        newChart.update('none');
      }
    }
  }

  function renderChartAngular(hist, theory) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) return;
    // 如当前不是柱状图，销毁后用柱状图重建
    if (state.chart && (!state.chart.config || state.chart.config.type !== 'bar')) {
      try { state.chart.destroy(); } catch (e) { console.warn('图表销毁失败:', e); }
      state.chart = null;
    }
    const bins = hist.counts.length;
    const centers = (theory && theory.centers) || (() => { const n = hist.counts.length; const a = new Array(n); for (let i = 0; i < n; i++) a[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]); return a; })();

    const chart = ensureChart();
    if (!chart) return;

    // 捕获角向理论曲线的可见性状态
    if (chart.data && chart.data.datasets) {
      const theoryIndex = chart.data.datasets.findIndex(ds => ds.label && ds.label.includes('理论曲线'));
      if (theoryIndex !== -1) {
        state.angularTheoryHidden = !chart.isDatasetVisible(theoryIndex);
      }
    }

    chart.data.labels = centers.map(v => v.toFixed(2));
    const datasets = [
      {
        label: '角向概率密度 (归一化)',
        data: Array.from(hist.counts),
        backgroundColor: 'rgba(255,255,255,0.7)',
        borderColor: 'rgba(255,255,255,0.95)',
        borderWidth: 1,
        borderRadius: 2,
        borderSkipped: false,
        order: 10
      }
    ];

    if (theory && theory.values && theory.values.length) {
      datasets.push({
        type: 'line',
        label: '理论曲线 (sin θ)',
        data: theory.values.map((y, i) => ({ x: i, y })),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        pointHoverRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y',
        tension: 0.2,
        order: 0,
        hidden: state.angularTheoryHidden
      });
    }
    applyThinLineStyle(datasets);
    chart.data.datasets = datasets;

    // 强制设置全局元素配置，确保线条不显示点
    if (!chart.options.elements) chart.options.elements = {};
    chart.options.elements.point = { radius: 0, hoverRadius: 0 };

    // 添加坐标轴标题和范围限制
    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = {
          display: true,
          text: '角度 θ (弧度)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
        // 优化x轴显示
        chart.options.scales.x.ticks.autoSkip = true;
        chart.options.scales.x.ticks.maxTicksLimit = 15;
        delete chart.options.scales.x.ticks.callback;
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = {
          display: true,
          text: '概率密度 P(θ)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
        // 【新增】固定Y轴范围 0-1
        chart.options.scales.y.min = 0;
        chart.options.scales.y.max = 1;
      }

      // 角向图不需要 y2，防止从能量图残留
      delete chart.options.scales.y2;
    }

    applyRollingAnimationPolicy(chart);
    // angular/phi 固定 y=0..1；无需 y 轴锁定
    chart.update('none'); // 使用无动画更新，避免曲线跳动
  }

  function renderChartPhi(hist, theory) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) return;
    // 如当前不是柱状图，销毁后用柱状图重建
    if (state.chart && (!state.chart.config || state.chart.config.type !== 'bar')) {
      try { state.chart.destroy(); } catch (e) { console.warn('图表销毁失败:', e); }
      state.chart = null;
    }
    const bins = hist.counts.length;
    const centers = (theory && theory.centers) || (() => { const n = hist.counts.length; const a = new Array(n); for (let i = 0; i < n; i++) a[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]); return a; })();

    const chart = ensureChart();
    if (!chart) return;

    // 捕获φ角向理论曲线的可见性状态
    if (chart.data && chart.data.datasets) {
      const theoryIndex = chart.data.datasets.findIndex(ds => ds.label && ds.label.includes('理论曲线'));
      if (theoryIndex !== -1) {
        state.phiTheoryHidden = !chart.isDatasetVisible(theoryIndex);
      }
    }

    chart.data.labels = centers.map(v => v.toFixed(2));
    const datasets = [
      {
        label: '方位角概率密度 (归一化)',
        data: Array.from(hist.counts),
        backgroundColor: 'rgba(255,255,255,0.7)',
        borderColor: 'rgba(255,255,255,0.95)',
        borderWidth: 1,
        borderRadius: 2,
        borderSkipped: false,
        order: 10
      }
    ];

    if (theory && theory.values && theory.values.length) {
      datasets.push({
        type: 'line',
        label: '理论曲线',
        data: theory.values.map((y, i) => ({ x: i, y })),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        pointHoverRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y',
        tension: 0.2,
        order: 0,
        hidden: state.phiTheoryHidden
      });
    }
    applyThinLineStyle(datasets);
    chart.data.datasets = datasets;

    // 【关键修复】φ角向图也不需要 y2，防止从能量图残留
    if (chart.options.scales) {
      delete chart.options.scales.y2;
    }

    // 强制设置全局元素配置，确保线条不显示点
    if (!chart.options.elements) chart.options.elements = {};
    chart.options.elements.point = { radius: 0, hoverRadius: 0 };

    // 添加坐标轴标题和范围限制
    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = {
          display: true,
          text: '角度 φ (弧度)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
        // 优化x轴显示
        chart.options.scales.x.ticks.autoSkip = true;
        chart.options.scales.x.ticks.maxTicksLimit = 15;
        delete chart.options.scales.x.ticks.callback;
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = {
          display: true,
          text: '概率密度 P(φ)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
        // 【新增】固定Y轴范围 0-1
        chart.options.scales.y.min = 0;
        chart.options.scales.y.max = 1;
      }
    }

    applyRollingAnimationPolicy(chart);
    // angular/phi 固定 y=0..1；无需 y 轴锁定
    chart.update('none'); // 使用无动画更新，避免曲线跳动
  }

  function renderChartPotential(theory) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) {
      return;
    }

    // 保存理论曲线可见性状态
    if (state.chart && state.chart.data && state.chart.data.datasets) {
      const theoryIndex = state.chart.data.datasets.findIndex(ds => ds.label === '理论积分曲线');
      if (theoryIndex !== -1) {
        state.potentialTheoryHidden = !state.chart.isDatasetVisible(theoryIndex);
      }
    }

    // 若当前不是柱状图，需销毁后重新创建
    if (state.chart && state.chart.config && state.chart.config.type !== 'bar') {
      try { state.chart.destroy(); } catch (e) { console.warn('图表销毁失败:', e); }
      state.chart = null;
    }

    const chart = ensureChart();
    if (!chart) {
      return;
    }

    // 准备 X 轴标签和数据集
    const datasets = [];
    let labels = [];

    // 理论曲线 - 折线叠加
    if (theory && theory.points && theory.points.length > 0) {
      labels = theory.points.map(p => p.x.toFixed(2));

      datasets.push({
        type: 'line',
        label: '理论积分曲线',
        data: theory.points.map((p, i) => ({ x: i, y: p.y })),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y',
        tension: 0.2,
        order: 5,
        hidden: state.potentialTheoryHidden || false
      });
    }



    applyThinLineStyle(datasets);
    chart.data.labels = labels;
    chart.data.datasets = datasets;

    // 更新坐标轴标题
    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = {
          display: true,
          text: '距离 r (a₀)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
        chart.options.scales.x.ticks.autoSkip = true;
        chart.options.scales.x.ticks.maxTicksLimit = 15;
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = {
          display: true,
          text: '累积势能 V(r) (Hartree)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
      }

    }

    applyRollingAnimationPolicy(chart);
    applyRollingYAxisLock(chart, datasets, 'bar:potential', 'potential');
    chart.update('none');
  }

  // 渲染势能密度 dE/dr 曲线（类似 renderChartPotential）
  function renderChartDEdr(theory) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) {
      return;
    }

    // 保存理论曲线可见性状态
    if (state.chart && state.chart.data && state.chart.data.datasets) {
      const theoryIndex = state.chart.data.datasets.findIndex(ds => ds.label === '理论 dE/dr');
      if (theoryIndex !== -1) {
        state.dEdrTheoryHidden = !state.chart.isDatasetVisible(theoryIndex);
      }
    }

    // 若当前不是柱状图，需销毁后重新创建
    if (state.chart && state.chart.config && state.chart.config.type !== 'bar') {
      try { state.chart.destroy(); } catch (e) { console.warn('图表销毁失败:', e); }
      state.chart = null;
    }

    const chart = ensureChart();
    if (!chart) {
      return;
    }

    // 准备 X 轴标签和数据集
    const datasets = [];
    let labels = [];

    // 理论曲线 - 折线叠加
    if (theory && theory.points && theory.points.length > 0) {
      labels = theory.points.map(p => p.x.toFixed(2));

      datasets.push({
        type: 'line',
        label: '理论 dE/dr',
        data: theory.points.map((p, i) => ({ x: i, y: p.y })),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y',
        tension: 0.2,
        order: 5,
        hidden: state.dEdrTheoryHidden || false
      });
    }

    // Z_eff 曲线 - 绿色，右侧坐标轴
    if (theory && theory.zeff && theory.zeff.length) {
      datasets.push({
        type: 'line',
        label: '有效核电荷 Z_eff',
        data: theory.zeff.map((y, i) => ({ x: i, y })),
        borderColor: 'rgba(0, 220, 0, 0.9)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y2',
        borderDash: [2, 2],
        tension: 0.2,
        order: 0,
        hidden: state.zeffHidden || false
      });
    }

    applyThinLineStyle(datasets);
    chart.data.labels = labels;
    chart.data.datasets = datasets;

    // 更新坐标轴标题
    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = {
          display: true,
          text: '距离 r (a₀)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
        chart.options.scales.x.ticks.autoSkip = true;
        chart.options.scales.x.ticks.maxTicksLimit = 15;
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = {
          display: true,
          text: '势能密度 dV/dr (Hartree/a₀)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
      }

      // Z_eff 坐标轴 - 右侧绿色
      if (theory && theory.zeff && theory.zeff.length) {
        const maxZ = Math.ceil(Math.max(...theory.zeff, 1));
        chart.options.scales.y2 = {
          position: 'right',
          grid: { drawOnChartArea: false },
          min: 0,
          max: maxZ > 1 ? maxZ : 2,
          ticks: {
            display: true,
            color: 'rgba(0, 220, 0, 0.9)',
            font: { size: 9 }
          },
          title: {
            display: true,
            text: 'Z_eff',
            color: 'rgba(0, 220, 0, 0.9)',
            font: { size: 10 }
          }
        };
      } else {
        // 没有 Z_eff 时移除 y2，避免残留轴
        delete chart.options.scales.y2;
      }
    }

    applyRollingAnimationPolicy(chart);
    applyRollingYAxisLock(chart, datasets, 'bar:dEdr', 'dEdr');
    chart.update('none');
  }

  // 渲染能量期望密度图 ε·P(r) - 样式照抄径向分布
  function renderChartEnergyDensity(hist, theory, options = {}) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) return;

    // 销毁非柱状图
    if (state.chart && state.chart.config && state.chart.config.type !== 'bar') {
      try { state.chart.destroy(); } catch (e) { }
      state.chart = null;
    }

    const chart = ensureChart();
    if (!chart) return;

    // 直方图中心作为 x 轴
    const centers = (theory && theory.centers) || (() => {
      const n = hist.counts.length;
      const a = new Array(n);
      for (let i = 0; i < n; i++) a[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);
      return a;
    })();

    // 动态调整显示精度
    const maxValue = Math.max(...centers);
    const decimalPlaces = maxValue > 100 ? 1 : (maxValue > 10 ? 2 : 3);

    chart.data.labels = centers.map(v => v.toFixed(decimalPlaces));
    const datasets = [
      {
        label: '能量期望密度 ε·P(r) (采样)',
        data: Array.from(hist.counts),
        backgroundColor: 'rgba(255,255,255,0.7)',
        borderColor: 'rgba(255,255,255,0.95)',
        borderWidth: 0,
        barPercentage: 1.0,
        categoryPercentage: 1.0,
        borderRadius: 2,
        borderSkipped: false,
        order: 10,
      },
    ];

    // 理论曲线
    if (theory && theory.values && theory.values.length) {
      datasets.push({
        type: 'line',
        label: '理论曲线 ε·P(r)',
        data: theory.values.map((y, i) => ({ x: i, y })),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y',
        tension: 0.2,
        order: 5,
      });
    }

    applyThinLineStyle(datasets);
    chart.data.datasets = datasets;

    const lockKey = (options && options.lockKey) ? options.lockKey : 'bar:energyDensity';
    const typeHint = (options && options.typeHint) ? options.typeHint : 'energyDensity';
    const yTitleText = (options && options.yTitle) ? options.yTitle : '能量密度 ε·P(r) (Hartree/a₀)';

    // 坐标轴配置
    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = {
          display: true,
          text: '径向距离 r (a₀)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = {
          display: true,
          text: yTitleText,
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };

        // 【关键修复】从角向/φ图切换过来时，清除固定的 0-1 轴范围
        delete chart.options.scales.y.min;
        delete chart.options.scales.y.max;
      }

      // 能量密度图不需要 y2
      delete chart.options.scales.y2;
    }

    applyRollingAnimationPolicy(chart);
    applyRollingYAxisLock(chart, datasets, lockKey, typeHint);
    chart.update('none');
  }

  // 渲染能量期望密度图（gamma变换）：数据与 ε·P(r) 完全同源，只在绘图阶段对 y 做 gamma=0.1
  function renderChartEnergyDensityGamma(hist, theory, gamma = 0.1) {
    const safeGamma = (typeof gamma === 'number' && Number.isFinite(gamma) && gamma > 0) ? gamma : 0.1;
    const gammaCorrect = (y) => {
      if (typeof y !== 'number' || !Number.isFinite(y)) return 0;
      if (y === 0) return 0;
      const s = y < 0 ? -1 : 1;
      const a = Math.abs(y);
      // 幂次 gamma：保符号，仅对幅值做变换
      const v = Math.pow(a, safeGamma);
      return Number.isFinite(v) ? (s * v) : 0;
    };

    const transformedHist = {
      ...hist,
      counts: Array.isArray(hist?.counts) ? hist.counts.map(gammaCorrect) : Array.from(hist?.counts || [], gammaCorrect)
    };
    const transformedTheory = {
      ...theory,
      values: Array.isArray(theory?.values) ? theory.values.map(gammaCorrect) : Array.from(theory?.values || [], gammaCorrect)
    };

    renderChartEnergyDensity(transformedHist, transformedTheory, {
      lockKey: 'bar:energyDensityGamma',
      typeHint: 'energyDensityGamma',
      yTitle: '能量密度 ε·P(r)（γ=0.1）'
    });
  }

  // 渲染能量期望累计图 E(R) = ∫ε·P(r)dr - 样式照抄径向分布
  // 诊断：累计能量曲线出现“大阶跃/长平台”时，判断是否属于“理论上本就期望命中极少（空 bin 合理）”，
  // 还是“理论期望命中数足够大但采样为 0（疑似采样/截断/链路错误）”。
  // 注意：该诊断必须基于未做 gamma 的累计能量 E(R)；gamma 后的 E(R)^γ 不再与概率质量线性对应。
  function diagEnergyCumulativeSummary(args) {
    try {
      const tag = args && args.tag;
      const lockKeyForDiag = args && args.lockKeyForDiag;
      const logPrefix = (args && args.logPrefix) ? String(args.logPrefix) : '';
      const cum = Array.from((args && args.cum) || []);
      const theoryCum = Array.from((args && args.theoryCum) || []);
      const edges = (args && args.edges) || [];
      const sampleCount = (args && Number.isFinite(args.sampleCount)) ? args.sampleCount : (
        (window.ElectronCloud && window.ElectronCloud.state && Array.isArray(window.ElectronCloud.state.radialSamples))
          ? window.ElectronCloud.state.radialSamples.length
          : 0
      );

      if (cum.length < 2) return null;
      const n = theoryCum.length ? Math.min(cum.length, theoryCum.length) : cum.length;
      if (n < 2) return null;

      if (!state._energyCumulativeDiagSeen) state._energyCumulativeDiagSeen = {};
      const seenKey = `${lockKeyForDiag || ''}|${tag || ''}|${n}`;
      if (state._energyCumulativeDiagSeen[seenKey]) return null;
      state._energyCumulativeDiagSeen[seenKey] = true;

      const cumEnd = cum[n - 1];
      const scale = Math.max(1e-12, Math.abs(cumEnd));
      const tolFlat = Math.max(1e-14, scale * 1e-8);

      let flatBins = 0;
      let longestFlat = 0;
      let currentFlat = 0;
      let longestFlatEnd = -1;

      // theoryCum[end] 应该趋近 ε（符号可能为负）
      const epsilonRef = (theoryCum.length && Number.isFinite(theoryCum[n - 1]) && Math.abs(theoryCum[n - 1]) > 1e-14)
        ? theoryCum[n - 1]
        : null;

      // 采样侧的 ε 参考（用于把采样累计增量换算回概率质量）
      const epsilonObsRef = (Number.isFinite(cum[n - 1]) && Math.abs(cum[n - 1]) > 1e-14)
        ? cum[n - 1]
        : null;

      let suspiciousExpectedBins = 0;
      let longestSuspiciousExpected = 0;
      let currentSuspiciousExpected = 0;
      let longestSuspiciousExpectedEnd = -1;
      let maxExpectedCount = 0;

      for (let i = 1; i < n; i++) {
        const d = cum[i] - cum[i - 1];
        const isFlat = Math.abs(d) <= tolFlat;
        if (isFlat) {
          flatBins++;
          currentFlat++;
          if (currentFlat > longestFlat) {
            longestFlat = currentFlat;
            longestFlatEnd = i;
          }
        } else {
          currentFlat = 0;
        }

        if (theoryCum.length && epsilonRef !== null && Number.isFinite(sampleCount) && sampleCount > 0) {
          const dt = theoryCum[i] - theoryCum[i - 1];
          // probMass = |ΔE / ε|
          const probMass = Math.abs(dt / epsilonRef);
          const expected = sampleCount * probMass;
          if (Number.isFinite(expected)) {
            if (expected > maxExpectedCount) maxExpectedCount = expected;
            const isSuspiciousExpected = isFlat && expected >= 3;
            if (isSuspiciousExpected) {
              suspiciousExpectedBins++;
              currentSuspiciousExpected++;
              if (currentSuspiciousExpected > longestSuspiciousExpected) {
                longestSuspiciousExpected = currentSuspiciousExpected;
                longestSuspiciousExpectedEnd = i;
              }
            } else {
              currentSuspiciousExpected = 0;
            }
          }
        }
      }

      // 近核专项统计：用户关心的“近核空缺”主要发生在 r 很小的若干 bin。
      // 这里不武断判错，只给出：在 r<=nearCoreMaxR 区间内，理论期望命中数/采样估计命中数与“高期望但实际为 0”的最坏情况。
      const nearCoreMaxR = (args && Number.isFinite(args.nearCoreMaxR)) ? args.nearCoreMaxR : 2.0;
      let nearCoreBins = 0;
      let nearCoreExpectedTotal = 0;
      let nearCoreExpectedMax = 0;
      let nearCoreObsApproxTotal = 0;
      let nearCoreZeroBinsExpectedGE0_5 = 0;
      let nearCoreWorstZeroExpected = 0;
      let nearCoreWorstZeroR = null;

      if (edges && edges.length >= n + 1 && epsilonRef !== null && Number.isFinite(sampleCount) && sampleCount > 0) {
        for (let i = 1; i < n; i++) {
          const rRight = edges[i + 0];
          if (!Number.isFinite(rRight) || rRight > nearCoreMaxR) break;
          nearCoreBins++;

          const dt = theoryCum.length ? (theoryCum[i] - theoryCum[i - 1]) : 0;
          const probMassTh = Math.abs(dt / epsilonRef);
          const expected = sampleCount * probMassTh;
          if (Number.isFinite(expected)) {
            nearCoreExpectedTotal += expected;
            if (expected > nearCoreExpectedMax) nearCoreExpectedMax = expected;
          }

          let obsApprox = null;
          if (epsilonObsRef !== null) {
            const dObs = cum[i] - cum[i - 1];
            const probMassObs = Math.abs(dObs / epsilonObsRef);
            obsApprox = sampleCount * probMassObs;
            if (Number.isFinite(obsApprox)) nearCoreObsApproxTotal += obsApprox;

            const isObsZero = Math.abs(dObs) <= tolFlat;
            if (isObsZero && Number.isFinite(expected) && expected >= 0.5) {
              nearCoreZeroBinsExpectedGE0_5++;
              if (expected > nearCoreWorstZeroExpected) {
                nearCoreWorstZeroExpected = expected;
                nearCoreWorstZeroR = [edges[i - 1], edges[i]];
              }
            }
          }
        }
      }

      const flatStart = (longestFlatEnd >= 0) ? (longestFlatEnd - longestFlat + 1) : -1;
      const suspiciousExpectedStart = (longestSuspiciousExpectedEnd >= 0)
        ? (longestSuspiciousExpectedEnd - longestSuspiciousExpected + 1)
        : -1;
      const rRange = (start, end) => {
        if (!edges || edges.length < end + 2 || start < 0 || end < 0) return null;
        return [edges[start], edges[end + 1]];
      };

      const payload = {
        tag,
        bins: n,
        cumEnd,
        sampleCount,
        flatBins,
        flatRatio: Number.isFinite(flatBins / (n - 1)) ? Number((flatBins / (n - 1)).toFixed(3)) : null,
        longestFlat,
        longestFlatR: rRange(flatStart, longestFlatEnd),
        suspiciousExpectedBins,
        longestSuspiciousExpected,
        longestSuspiciousExpectedR: rRange(suspiciousExpectedStart, longestSuspiciousExpectedEnd),
        maxExpectedCount: Number.isFinite(maxExpectedCount) ? Number(maxExpectedCount.toFixed(2)) : null,

        nearCoreMaxR,
        nearCoreBins,
        nearCoreExpectedTotal: Number.isFinite(nearCoreExpectedTotal) ? Number(nearCoreExpectedTotal.toFixed(3)) : null,
        nearCoreExpectedMax: Number.isFinite(nearCoreExpectedMax) ? Number(nearCoreExpectedMax.toFixed(3)) : null,
        nearCoreObsApproxTotal: Number.isFinite(nearCoreObsApproxTotal) ? Number(nearCoreObsApproxTotal.toFixed(3)) : null,
        nearCoreZeroBinsExpectedGE0_5,
        nearCoreWorstZeroExpected: Number.isFinite(nearCoreWorstZeroExpected) ? Number(nearCoreWorstZeroExpected.toFixed(3)) : null,
        nearCoreWorstZeroR
      };

      window.__EC_LAST_ENERGY_CUMULATIVE_DIAG = payload;

      const json = (() => {
        try { return JSON.stringify(payload); } catch (e) { return null; }
      })();

      const suffix = logPrefix ? ` ${logPrefix}` : '';
      if (suspiciousExpectedBins > 0) {
        console.warn(`[EnergyCumulativeDiag${suffix}]`, payload);
        if (json) console.warn(`[EnergyCumulativeDiag${suffix} JSON]`, json);
      } else {
        console.log(`[EnergyCumulativeDiag${suffix}]`, payload);
        if (json) console.log(`[EnergyCumulativeDiag${suffix} JSON]`, json);
      }

      return payload;
    } catch (e) {
      return null;
    }
  }

  function renderChartEnergyCumulative(hist, theory, options = {}) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) return;

    // 销毁非柱状图
    if (state.chart && state.chart.config && state.chart.config.type !== 'bar') {
      try { state.chart.destroy(); } catch (e) { }
      state.chart = null;
    }

    const chart = ensureChart();
    if (!chart) return;

    // 直方图中心作为 x 轴
    const centers = (theory && theory.centers) || (() => {
      const n = hist.counts.length;
      const a = new Array(n);
      for (let i = 0; i < n; i++) a[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);
      return a;
    })();

    // 动态调整显示精度
    const maxValue = Math.max(...centers);
    const decimalPlaces = maxValue > 100 ? 1 : (maxValue > 10 ? 2 : 3);

    const lockKey = (options && options.lockKey) ? options.lockKey : 'bar:energyCumulative';
    const typeHint = (options && options.typeHint) ? options.typeHint : 'energyCumulative';
    const yTitleText = (options && options.yTitle) ? options.yTitle : '累计能量 E(R) (Hartree)';
    const samplingLabel = (options && options.samplingLabel) ? options.samplingLabel : '能量期望累计 E(R) (采样)';
    const theoryLabel = (options && options.theoryLabel) ? options.theoryLabel : '理论曲线 E(R) → ε';

    // 仅对累计能量图做诊断（可通过 options.skipDiag 禁用）
    if (!(options && options.skipDiag)) {
      diagEnergyCumulativeSummary({
        tag: typeHint,
        lockKeyForDiag: lockKey,
        cum: Array.from(hist?.counts || []),
        edges: hist?.edges || [],
        theoryCum: Array.from(theory?.values || [])
      });
    }

    chart.data.labels = centers.map(v => v.toFixed(decimalPlaces));
    const datasets = [
      {
        label: samplingLabel,
        data: Array.from(hist.counts),
        backgroundColor: 'rgba(255,255,255,0.7)',
        borderColor: 'rgba(255,255,255,0.95)',
        borderWidth: 0,
        barPercentage: 1.0,
        categoryPercentage: 1.0,
        borderRadius: 2,
        borderSkipped: false,
        order: 10,
      },
    ];

    // 理论曲线
    if (theory && theory.values && theory.values.length) {
      datasets.push({
        type: 'line',
        label: theoryLabel,
        data: theory.values.map((y, i) => ({ x: i, y })),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        borderWidth: 2.5,
        yAxisID: 'y',
        tension: 0.2,
        order: 5,
      });
    }

    applyThinLineStyle(datasets);
    chart.data.datasets = datasets;

    // 坐标轴配置
    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = {
          display: true,
          text: '径向距离 r (a₀)',
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = {
          display: true,
          text: yTitleText,
          color: '#d0d0d0',
          font: { size: 12, weight: '500' }
        };

        // 【关键修复】从角向/φ图切换过来时，清除固定的 0-1 轴范围
        delete chart.options.scales.y.min;
        delete chart.options.scales.y.max;
      }

      // 累计能量图不需要 y2
      delete chart.options.scales.y2;
    }

    applyRollingAnimationPolicy(chart);
    applyRollingYAxisLock(chart, datasets, lockKey, typeHint);
    chart.update('none');
  }

  // 渲染累计能量 E(R)（gamma变换）：数据与 E(R) 完全同源，只在绘图阶段对 y 做 gamma=0.1
  function renderChartEnergyCumulativeGamma(hist, theory, gamma = 0.1) {
    const safeGamma = (typeof gamma === 'number' && Number.isFinite(gamma) && gamma > 0) ? gamma : 0.1;
    const gammaCorrect = (y) => {
      if (typeof y !== 'number' || !Number.isFinite(y)) return 0;
      if (y === 0) return 0;
      const s = y < 0 ? -1 : 1;
      const a = Math.abs(y);
      const v = Math.pow(a, safeGamma);
      return Number.isFinite(v) ? (s * v) : 0;
    };

    // 诊断：gamma 图表必须基于 preGamma 的累计 E(R) 来做，否则“理论期望命中数”不再成立
    // 为了不污染 UI，这里仅打印诊断并保存全局变量。
    try {
      const payload = diagEnergyCumulativeSummary({
        tag: 'energyCumulativeGamma-preGamma',
        lockKeyForDiag: 'bar:energyCumulativeGamma:preGamma',
        logPrefix: 'preGamma',
        cum: Array.from(hist?.counts || []),
        edges: hist?.edges || [],
        theoryCum: Array.from(theory?.values || [])
      });
      if (payload) window.__EC_LAST_ENERGY_CUMULATIVE_PRE_GAMMA_DIAG = payload;
    } catch (e) { }

    const transformedHist = {
      ...hist,
      counts: Array.isArray(hist?.counts) ? hist.counts.map(gammaCorrect) : Array.from(hist?.counts || [], gammaCorrect)
    };
    const transformedTheory = {
      ...theory,
      values: Array.isArray(theory?.values) ? theory.values.map(gammaCorrect) : Array.from(theory?.values || [], gammaCorrect)
    };

    renderChartEnergyCumulative(transformedHist, transformedTheory, {
      lockKey: 'bar:energyCumulativeGamma',
      typeHint: 'energyCumulativeGamma',
      yTitle: '累计能量 E(R)（γ=0.1）',
      samplingLabel: '累计能量 E(R)（γ=0.1，采样）',
      theoryLabel: '理论曲线 E(R)（γ=0.1）',
      // 禁用对 gamma 后曲线的诊断（gamma 后不再与概率质量线性对应）
      skipDiag: true
    });
  }

  // 保留旧函数名以兼容
  function renderChartEnergy(theory) {
    // 兼容旧调用：使用 energyDensity
    renderChartEnergyDensity({ counts: theory.dEdr || [], edges: [] }, theory);
  }
  function renderChartLocalEnergy(theory) {
    renderChartEnergy(theory);
  }


  // 渲染能量密度 dE/dr 双对数曲线
  function renderChartDEdrLog(theory) {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) return;

    // 保存理论曲线可见性状态
    if (state.chart && state.chart.data && state.chart.data.datasets) {
      const theoryIndex = state.chart.data.datasets.findIndex(ds => ds.label === '理论 dE/dr');
      if (theoryIndex !== -1) {
        state.dEdrLogTheoryHidden = !state.chart.isDatasetVisible(theoryIndex);
      }
    }

    // 使用折线图
    if (!ensureChartLine()) return;
    const chart = state.chart;

    const datasets = [];
    let labels = [];

    if (theory && theory.points && theory.points.length > 0) {
      labels = theory.points.map(p => p.x.toFixed(2));
      datasets.push({
        type: 'line',
        label: '理论 dE/dr',
        // 【关键修复】与直方图对齐
        data: theory.points.map(p => p.y),
        borderColor: 'rgba(255, 255, 255, 0.95)',
        backgroundColor: 'transparent',
        pointRadius: 0,
        borderWidth: 2.0,
        tension: 0.4, // 平滑
        cubicInterpolationMode: 'monotone',
        order: 5,
        hidden: state.dEdrLogTheoryHidden || false
      });
    }

    applyThinLineStyle(datasets);
    chart.data.labels = labels;
    chart.data.datasets = datasets;

    if (chart.options.scales) {
      if (chart.options.scales.x) {
        chart.options.scales.x.title = { display: true, text: 'log₁₀(r) (r in a₀)', color: '#d0d0d0', font: { size: 12 } };
      }
      if (chart.options.scales.y) {
        chart.options.scales.y.title = { display: true, text: '势能密度 dV/dr (Hartree/a₀) [对数刻度]', color: '#d0d0d0', font: { size: 12 } };
      }
    }

    chart.update('none');
  }

  // 确保创建一个折线图实例
  function ensureChartLine() {
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart) return false;

    if (state.chart && state.chart.config.type === 'line') return true;
    if (state.chart) {
      state.chart.destroy();
      state.chart = null;
    }

    // 创建新图表
    state.chart = new Chart(ctx, {
      type: 'line',
      data: { datasets: [] },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: { duration: 0 },
        interaction: { mode: 'index', intersect: false },
        plugins: {
          legend: { display: true },
          tooltip: {
            enabled: true,
            callbacks: {
              label: function (context) {
                return context.dataset.label + ': ' + context.parsed.y.toFixed(4);
              }
            }
          }
        },
        scales: {
          x: {
            type: 'linear',
            position: 'bottom',
            title: { display: true, text: 'r (a0)', color: '#d0d0d0' },
            min: 0,
            ticks: { color: '#d0d0d0', maxTicksLimit: 10 },
            grid: { color: 'rgba(255,255,255,0.08)' }
          },
          y: {
            position: 'left',
            title: { display: true, text: 'Value', color: '#d0d0d0' },
            ticks: { color: '#d0d0d0' },
            grid: { color: 'rgba(255,255,255,0.08)' }
          }
        }
      }
    });
    return true;
  }

  // 对比模式专用：渲染散点图
  function renderChartCompare(orbitalDataMap, type) {
    if (typeof type === 'string') type = type.trim(); // 自动清理空格
    console.log('renderChartCompare 被调用，类型:', type, 'orbitalDataMap:', orbitalDataMap ? Object.keys(orbitalDataMap).length : 0); // 调试信息
    const ctx = document.getElementById('probability-chart');
    if (!ctx || !window.Chart || !window.Hydrogen) {
      console.log('图表渲染失败：canvas、Chart.js或Hydrogen不可用'); // 调试信息
      return;
    }

    // 【修复】在任何图表存在时都尝试捕获理论曲线状态
    // 检查当前图表中是否有理论曲线数据集，保存其可见性状态
    if (state.chart && state.chart.data && state.chart.data.datasets) {
      const theoryIndex = state.chart.data.datasets.findIndex(ds => ds.label && ds.label.includes('理论'));
      if (theoryIndex !== -1) {
        // isDatasetVisible返回true表示可见，取反得到hidden状态
        state.compareTheoryHidden = !state.chart.isDatasetVisible(theoryIndex);
      }
    }

    const rollingEnabled = syncRollingCache();

    // 如果当前不是折线图，销毁后重新创建
    if (state.chart && state.chart.config && state.chart.config.type !== 'line') {
      console.log('销毁旧图表，创建折线图'); // 调试信息
      try { state.chart.destroy(); } catch (e) { console.warn('图表销毁失败:', e); }
      state.chart = null;
    }

    // 创建折线图
    if (!state.chart) {
      console.log('创建新的折线图实例'); // 调试信息
      ctx.style.backgroundColor = 'transparent';

      state.chart = new Chart(ctx, {
        type: 'line',
        data: { datasets: [] },
        options: {
          responsive: true,
          maintainAspectRatio: false,
          backgroundColor: 'transparent',
          animation: {
            duration: rollingEnabled ? 0 : 300,
            easing: 'easeOutQuart'
          },
          interaction: {
            intersect: false,
            mode: 'index'
          },
          layout: {
            padding: 0
          },
          plugins: {
            legend: {
              labels: {
                color: '#e8e8e8',
                font: { size: 11, weight: '500' },
                padding: 12
              }
            },
            tooltip: {
              backgroundColor: 'rgba(20, 20, 20, 0.95)',
              titleColor: '#ffffff',
              bodyColor: '#e8e8e8',
              borderColor: 'rgba(255, 255, 255, 0.2)',
              borderWidth: 1,
              cornerRadius: 6,
              callbacks: {
                label: function (context) {
                  return `${context.dataset.label}: (${context.parsed.x.toFixed(3)}, ${context.parsed.y.toFixed(6)})`;
                }
              }
            }
          },
          scales: {
            x: {
              type: 'linear',
              min: (type === 'angular') ? 0 : ((type === 'phi' || type === 'azimuthal') ? 0 : undefined),
              max: (type === 'angular') ? Math.PI : ((type === 'phi' || type === 'azimuthal') ? (2 * Math.PI) : undefined),
              title: {
                display: true,
                text: (type === 'radial' || type === 'energyDensity' || type === 'energyDensityGamma' || type === 'energyCumulative' || type === 'energyCumulativeGamma') ? '径向距离 r (a₀)' :
                  (type === 'angular') ? '角度 θ (弧度)' :
                    (type === 'phi' || type === 'azimuthal') ? '角度 φ (弧度)' :
                      (type === 'potential' || type === 'dEdr' || type === 'localEnergy') ? '距离 r (a₀)' : '角度 (弧度)',
                color: '#d0d0d0',
                font: { size: 12, weight: '500' }
              },
              ticks: {
                color: '#d0d0d0',
                font: { size: 10 }
              },
              grid: {
                color: 'rgba(255,255,255,0.08)',
                lineWidth: 1
              },
              border: {
                color: 'rgba(255,255,255,0.15)'
              }
            },
            y: {
              type: 'linear',
              title: {
                display: true,
                text: type === 'potential' ? '累积势能 V(r) (Hartree)' :
                  type === 'dEdr' ? '势能密度 dV/dr (Hartree/a₀)' :
                    type === 'localEnergy' ? '径向能量密度 (Hartree/a₀)' :
                      type === 'energyDensity' ? '能量密度 ε·P(r) (Hartree/a₀)' :
                        type === 'energyDensityGamma' ? '能量密度 ε·P(r)（γ=0.1）' :
                        type === 'energyCumulative' ? '累计能量 E(R) (Hartree)' :
                        type === 'energyCumulativeGamma' ? '累计能量 E(R)（γ=0.1）' :
                          '概率密度',
                color: '#d0d0d0',
                font: { size: 12, weight: '500' }
              },
              // 【新增】角向分布图固定Y轴范围 0-1
              min: (type === 'angular' || type === 'phi' || type === 'azimuthal') ? 0 : undefined,
              max: (type === 'angular' || type === 'phi' || type === 'azimuthal') ? 1 : undefined,
              ticks: {
                color: '#d0d0d0',
                font: { size: 10 }
              },
              grid: {
                color: 'rgba(255,255,255,0.08)',
                lineWidth: 1
              },
              border: {
                color: 'rgba(255,255,255,0.15)'
              }
            }
          },
          elements: {
            point: {
              radius: 0, // 不显示数据点
              hoverRadius: 3
            },
            line: {
              tension: 0.1 // 轻微的曲线平滑
            }
          }
        }
      });
    } else {
      // 如果图表已存在，更新坐标轴标题
      if (state.chart.options.scales) {
        if (state.chart.options.scales.x && state.chart.options.scales.x.title) {
          if (type === 'potentialLog' || type === 'dEdrLog') {
            // 已移除 log-log 图表类型
          } else if (
            type === 'radial' ||
            type === 'potential' ||
            type === 'dEdr' ||
            type === 'localEnergy' ||
            type === 'energyDensity' ||
            type === 'energyDensityGamma' ||
            type === 'energyCumulative' ||
            type === 'energyCumulativeGamma'
          ) {
            state.chart.options.scales.x.title.text = '距离 r (a₀)';
          } else {
            state.chart.options.scales.x.title.text = '角度 (弧度)';
          }

          // 【关键修复】钉死角向图 x 轴范围
          if (type === 'angular') {
            state.chart.options.scales.x.min = 0;
            state.chart.options.scales.x.max = Math.PI;
          } else if (type === 'phi' || type === 'azimuthal') {
            state.chart.options.scales.x.min = 0;
            state.chart.options.scales.x.max = 2 * Math.PI;
          } else {
            delete state.chart.options.scales.x.min;
            delete state.chart.options.scales.x.max;
          }
        }
        if (state.chart.options.scales.y && state.chart.options.scales.y.title) {
          if (type === 'potential') {
            state.chart.options.scales.y.title.text = '累积势能 V(r) (Hartree)';
          } else if (type === 'dEdr') {
            state.chart.options.scales.y.title.text = '势能密度 dV/dr (Hartree/a₀)';
          } else if (type === 'localEnergy') {
            state.chart.options.scales.y.title.text = '径向能量密度 (Hartree/a₀)';
          } else if (type === 'energyDensity') {
            state.chart.options.scales.y.title.text = '能量密度 ε·P(r) (Hartree/a₀)';
          } else if (type === 'energyDensityGamma') {
            state.chart.options.scales.y.title.text = '能量密度 ε·P(r)（γ=0.1）';
          } else if (type === 'energyCumulative') {
            state.chart.options.scales.y.title.text = '累计能量 E(R) (Hartree)';
          } else if (type === 'energyCumulativeGamma') {
            state.chart.options.scales.y.title.text = '累计能量 E(R)（γ=0.1）';
          } else {
            state.chart.options.scales.y.title.text = '概率密度';
          }

          // 【关键修复】钉死角向图 y 轴范围
          if (type === 'angular' || type === 'phi' || type === 'azimuthal') {
            state.chart.options.scales.y.min = 0;
            state.chart.options.scales.y.max = 1;
          } else {
            // 非角向：通常让其自适应；但滚动模式会在后续统一锁定 y 轴
            delete state.chart.options.scales.y.min;
            delete state.chart.options.scales.y.max;
          }
        }
      }
    }

    // 准备数据集 - 先将原始数据转为直方图，再转为散点
    const datasets = [];

    // 使用统一的常量定义
    const compareColors = window.ElectronCloud?.constants?.compareColors || [
      { name: 'red', value: [1, 0.2, 0.2] },
      { name: 'green', value: [0.2, 1, 0.2] },
      { name: 'blue', value: [0.2, 0.2, 1] }
    ];

    // 使用统一的轨道显示名称映射
    const orbitalDisplayNameMap = window.ElectronCloud?.constants?.orbitalDisplayNames || {};

    // 为了保持一致性，使用与普通模式相同的参数计算直方图
    let maxDistance = 0;
    let totalSamples = 0;

    // 首先找到最大距离用于确定动态范围（优先采样数据；没有采样时回退到理论估算半径）
    for (const [orbitalKey, samples] of Object.entries(orbitalDataMap || {})) {
      if (!samples || samples.length === 0) continue;
      totalSamples += samples.length;
      if (
        type === 'radial' ||
        type === 'potential' ||
        type === 'dEdr' ||
        type === 'localEnergy' ||
        type === 'energyDensity' ||
        type === 'energyDensityGamma' ||
        type === 'energyCumulative' ||
        type === 'energyCumulativeGamma'
      ) {
        // 【性能修复】使用循环替代Math.max(...array)，避免大数组栈溢出
        for (let i = 0; i < samples.length; i++) {
          if (samples[i].r > maxDistance) {
            maxDistance = samples[i].r;
          }
        }
      }
    }

    // 获取当前选择的轨道顺序，确保颜色分配与选择顺序一致
    // 【重构】比照模式下使用activeSlots获取slot配置（包含原子类型）
    const activeSlots = window.ElectronCloud?.state?.compareMode?.activeSlots || [];
    const currentOrbitals = window.ElectronCloud?.state?.currentOrbitals || Object.keys(orbitalDataMap || {});
    console.log('当前轨道顺序:', currentOrbitals, 'activeSlots:', activeSlots); // 调试信息

    // 【关键修复】若尚未采样（maxDistance=0），仍要在比照模式显示理论曲线
    // 使用 estimateOrbitalRadius95 为每个 slot 估算一个可用的 r 范围
    if ((type === 'radial' || type === 'potential' || type === 'dEdr' || type === 'localEnergy' || type === 'energyDensity' || type === 'energyDensityGamma' || type === 'energyCumulative' || type === 'energyCumulativeGamma') && maxDistance <= 0) {
      const estimateFn = window.Hydrogen?.estimateOrbitalRadius95;
      if (estimateFn && activeSlots && activeSlots.length > 0) {
        for (const slot of activeSlots) {
          const atomType = slot.atom || 'H';
          const orbitalKey = slot.orbital || '1s';
          const r95 = estimateFn(atomType, orbitalKey);
          if (Number.isFinite(r95) && r95 > maxDistance) maxDistance = r95;
        }
      }
      // 仍无法估算时，回退到当前状态的采样边界/默认值
      if (maxDistance <= 0) {
        const st = window.ElectronCloud?.state;
        maxDistance = Math.max(15, st?.samplingBoundary || 0, st?.farthestDistance || 0);
      }
    }

    // 滚动模式：锁定径向网格（bins/rmax），避免每次滚动导致横轴网格跳变
    let lockedRGrid = null;
    if (type === 'radial' || type === 'potential' || type === 'dEdr' || type === 'localEnergy' || type === 'energyDensity' || type === 'energyDensityGamma' || type === 'energyCumulative' || type === 'energyCumulativeGamma') {
      const dynamicRmaxCandidate = Math.max(1, maxDistance * 1.08);
      const baseBins = 240;
      const sampleDensity = totalSamples / Math.max(1, dynamicRmaxCandidate);
      const adaptiveBinsCandidate = Math.min(400, Math.max(baseBins, Math.floor(sampleDensity * 0.5)));

      if (rollingEnabled) {
        lockedRGrid = state.rollingCache.compareRGridLocks[type];
        if (!lockedRGrid) {
          lockedRGrid = { rmax: dynamicRmaxCandidate, bins: adaptiveBinsCandidate };
          state.rollingCache.compareRGridLocks[type] = lockedRGrid;
        }
      } else {
        lockedRGrid = { rmax: dynamicRmaxCandidate, bins: adaptiveBinsCandidate };
      }
    }

    // 按照选择顺序处理轨道数据
    for (let colorIndex = 0; colorIndex < activeSlots.length; colorIndex++) {
      const slotConfig = activeSlots[colorIndex];
      const atomType = slotConfig.atom || 'H';
      // 构建与sampling.js中相同的键
      const sampleKey = `${atomType}_${slotConfig.orbital}_slot${slotConfig.slotIndex}`;
      const samples = orbitalDataMap[sampleKey];

      const hasSamples = !!(samples && samples.length > 0);

      // 【关键修复】颜色必须绑定到 UI slotIndex，而不是 activeSlots 的数组位置
      const slotColorIndex = Number.isInteger(slotConfig.slotIndex) ? slotConfig.slotIndex : colorIndex;
      const color = compareColors[(slotColorIndex % compareColors.length + compareColors.length) % compareColors.length];
      // 【关键修复】标签中显示原子类型
      const displayName = `${slotConfig.atom} ${orbitalDisplayNameMap[slotConfig.orbital] || slotConfig.orbital}`;
      let hist, centers;
      let potentialValues; // 存储势能值

      if (type === 'radial') {
        // 使用与普通模式相同的参数（滚动模式下固定 bins/rmax）
        const dynamicRmax = lockedRGrid ? lockedRGrid.rmax : Math.max(1, maxDistance * 1.08);
        const adaptiveBins = lockedRGrid ? lockedRGrid.bins : 240;

        if (hasSamples) {
          // 提取径向数据
          const radialData = samples.map(s => s.r);
          hist = window.Hydrogen.histogramRadialFromSamples(radialData, adaptiveBins, dynamicRmax, true, false);
        } else {
          // 无采样时：构造一个空直方图，仅用于生成 centers 与理论曲线
          const edges = new Array(adaptiveBins + 1);
          const counts = new Array(adaptiveBins).fill(0);
          for (let i = 0; i <= adaptiveBins; i++) edges[i] = (dynamicRmax * i) / adaptiveBins;
          hist = { edges, counts };
        }

        // 计算bin中心
        centers = new Array(adaptiveBins);
        for (let i = 0; i < adaptiveBins; i++) {
          centers[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);
        }
      } else if (type === 'angular') {
        // θ角向数据
        const angularBins = 180;
        // 强制使用固定范围 [0, π]
        const edges = new Array(angularBins + 1);
        const counts = new Array(angularBins).fill(0);
        for (let i = 0; i <= angularBins; i++) edges[i] = (Math.PI * i) / angularBins;

        if (hasSamples) {
          const angularData = samples.map(s => s.theta);
          // 手动统计
          const binWidth = Math.PI / angularBins;
          for (let val of angularData) {
            const idx = Math.floor(val / Math.PI * angularBins);
            if (idx >= 0 && idx < angularBins) counts[idx]++;
          }
          // 归一化
          const total = angularData.length || 1;
          for (let i = 0; i < angularBins; i++) counts[i] /= (total * binWidth);
        }
        hist = { edges, counts };

        // 计算bin中心
        centers = new Array(angularBins);
        for (let i = 0; i < angularBins; i++) {
          centers[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);
        }
      } else if (type === 'potential' || type === 'dEdr' || type === 'localEnergy' || type === 'energyDensity' || type === 'energyDensityGamma' || type === 'energyCumulative' || type === 'energyCumulativeGamma') {
        // 能量图：复用采样决定的 r 网格范围
        const dynamicRmax = lockedRGrid ? lockedRGrid.rmax : Math.max(1, maxDistance * 1.08);
        const adaptiveBins = lockedRGrid ? lockedRGrid.bins : 240;

        if (hasSamples) {
          const radialData = samples.map(s => s.r);
          hist = window.Hydrogen.histogramRadialFromSamples(radialData, adaptiveBins, dynamicRmax, true, false);
        } else {
          const edges = new Array(adaptiveBins + 1);
          const counts = new Array(adaptiveBins).fill(0);
          for (let i = 0; i <= adaptiveBins; i++) edges[i] = (dynamicRmax * i) / adaptiveBins;
          hist = { edges, counts };
        }

        centers = new Array(adaptiveBins);
        for (let i = 0; i < adaptiveBins; i++) {
          centers[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);
        }
      } else {
        // φ角向数据 (azimuthal)
        const phiBins = 180;
        // 强制使用固定范围 [0, 2π]（与单选模式 histogramPhiFromSamples 保持一致）
        const edges = new Array(phiBins + 1);
        const counts = new Array(phiBins).fill(0);
        for (let i = 0; i <= phiBins; i++) edges[i] = (2 * Math.PI * i) / phiBins;

        if (hasSamples) {
          const phiData = samples.map(s => s.phi);
          // 使用手动统计填充 counts，确保 range 一致
          for (let val of phiData) {
            // 归一化到 [0, 2π)
            while (val < 0) val += 2 * Math.PI;
            while (val >= 2 * Math.PI) val -= 2 * Math.PI;

            const idx = Math.floor(val / (2 * Math.PI) * phiBins);
            if (idx >= 0 && idx < phiBins) counts[idx]++;
          }
          // 归一化 density
          const binWidth = 2 * Math.PI / phiBins;
          const total = phiData.length || 1;
          for (let i = 0; i < phiBins; i++) counts[i] /= (total * binWidth);
        }
        hist = { edges, counts };

        // 计算bin中心
        centers = new Array(phiBins);
        for (let i = 0; i < phiBins; i++) {
          centers[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);
        }
      }

      // 将颜色值从[0,1]范围转换为[0,255]范围
      const colorValues = color.value.map(v => Math.round(v * 255));

      // 采样曲线：用于概率分布图和能量密度/累计图
      if (hasSamples && type !== 'potential' && type !== 'dEdr' && type !== 'localEnergy') {
        let samplingData;

        if (type === 'energyDensity' || type === 'energyDensityGamma' || type === 'energyCumulative' || type === 'energyCumulativeGamma') {
          // 能量图：采样数据 × ε
          const baseKey = (slotConfig.orbital || '').replace(/[xyz]/g, '').replace(/_.*/, '');
          const energies = window.SlaterBasis && window.SlaterBasis[atomType]
            ? window.SlaterBasis[atomType].energies
            : null;
          const epsilon = (energies && energies[baseKey] !== undefined) ? energies[baseKey] : -0.5;
          const dr = (hist.edges[1] - hist.edges[0]) || 0.01;

          const gammaCorrect = (y) => {
            if (typeof y !== 'number' || !Number.isFinite(y)) return 0;
            if (y === 0) return 0;
            const s = y < 0 ? -1 : 1;
            const v = Math.pow(Math.abs(y), 0.1);
            return Number.isFinite(v) ? (s * v) : 0;
          };

          if (type === 'energyDensity') {
            samplingData = centers.map((center, index) => ({
              x: center,
              y: hist.counts[index] * epsilon
            }));
          } else if (type === 'energyDensityGamma') {
            samplingData = centers.map((center, index) => ({
              x: center,
              y: gammaCorrect(hist.counts[index] * epsilon)
            }));
          } else if (type === 'energyCumulative') {
            // energyCumulative: 累积求和
            let cumSum = 0;
            samplingData = centers.map((center, index) => {
              cumSum += hist.counts[index] * epsilon * dr;
              return { x: center, y: cumSum };
            });
          } else if (type === 'energyCumulativeGamma') {
            // energyCumulativeGamma：与 energyCumulative 同源（先得到 E(R)，再对 y 做 gamma=0.1）
            let cumSum = 0;
            samplingData = centers.map((center, index) => {
              cumSum += hist.counts[index] * epsilon * dr;
              return { x: center, y: gammaCorrect(cumSum) };
            });
          } else {
            // 默认 radial 或 fallback
            // 【严格修正】只有明确是 radial 时才按 radial 处理，避免 fallback 掩盖错误
            if (type === 'radial') {
              samplingData = centers.map((center, index) => ({
                x: center,
                y: hist.counts[index]
              }));
            } else if (type === 'angular') {
              // θ角向采样直方图
              samplingData = centers.map((center, index) => ({
                x: center,
                y: hist.counts[index]
              }));
            } else if (type === 'phi' || type === 'azimuthal') {
              // φ方位角采样直方图
              samplingData = centers.map((center, index) => ({
                x: center,
                y: hist.counts[index]
              }));
            } else {
              console.error('renderChartCompare: 未知或未实现的图表类型用于采样曲线:', type);
              samplingData = [];
            }
          }
        } else {
          // 普通概率分布图 - 处理所有类型
          if (type === 'radial') {
            samplingData = centers.map((center, index) => ({
              x: center,
              y: hist.counts[index]
            }));
          } else if (type === 'angular') {
            samplingData = centers.map((center, index) => ({
              x: center,
              y: hist.counts[index]
            }));
          } else if (type === 'phi' || type === 'azimuthal') {
            samplingData = centers.map((center, index) => ({
              x: center,
              y: hist.counts[index]
            }));
          } else {
            // 理论上 potential/dEdr/localEnergy 已经被 excluded，不应该进这里
            // 但为了安全
            samplingData = [];
          }
        }

        samplingData.sort((a, b) => a.x - b.x);

        datasets.push({
          label: displayName,
          data: samplingData,
          borderColor: `rgba(${colorValues.join(',')}, 1.0)`,
          backgroundColor: `rgba(${colorValues.join(',')}, 0.1)`,
          borderWidth: 1.5,
          pointRadius: 0,
          pointHoverRadius: 3,
          fill: false,
          tension: 0.1,
        });
      }

      // 【新增】为该轨道添加理论曲线（虚线）
      const orbitalParams = window.Hydrogen?.orbitalParamsFromKey(slotConfig.orbital);
      if (orbitalParams) {
        // const atomType = slotConfig.atom || 'H'; // Moved to top
        const atomZ = window.SlaterBasis && window.SlaterBasis[atomType] ? window.SlaterBasis[atomType].Z : 1;
        let theoryData;

        if (type === 'radial') {
          // 径向理论曲线
          theoryData = centers.map(r => ({
            x: r,
            y: window.Hydrogen.radialPDF(orbitalParams.n, orbitalParams.l, r, atomZ, 1, atomType)
          }));
        } else if (type === 'angular') {
          // θ 角向理论曲线
          theoryData = centers.map(theta => ({
            x: theta,
            y: window.Hydrogen.angularPDF_Theta(orbitalParams.l, orbitalParams.angKey.m, theta)
          }));
        } else if (type === 'potential') {
          // 势能理论曲线
          // 使用 calculateCumulativeOrbitalEnergy 计算理论 V(r)
          const res = window.Hydrogen.calculateCumulativeOrbitalEnergy(orbitalParams.n, orbitalParams.l, atomZ, atomType, centers);
          theoryData = centers.map((r, i) => ({ x: r, y: (res && res.E && res.E.length > i) ? res.E[i] : 0 }));
        } else if (type === 'dEdr') {
          const res = window.Hydrogen.calculateCumulativeOrbitalEnergy(orbitalParams.n, orbitalParams.l, atomZ, atomType, centers);
          theoryData = centers.map((r, i) => ({
            x: r,
            y: (res && res.dEdr && res.dEdr.length > i) ? res.dEdr[i] : 0
          }));
        } else if (type === 'energyDensity') {
            // 能量密度理论曲线：与采样侧严格同定义：ε·P(r)
            const baseKey = (slotConfig.orbital || '').replace(/[xyz]/g, '').replace(/_.*/, '');
            const energies = window.SlaterBasis && window.SlaterBasis[atomType]
              ? window.SlaterBasis[atomType].energies
              : null;
            const epsilon = (energies && energies[baseKey] !== undefined) ? energies[baseKey] : -0.5;
            theoryData = centers.map(r => ({
              x: r,
              y: epsilon * window.Hydrogen.radialPDF(orbitalParams.n, orbitalParams.l, r, atomZ, 1, atomType)
            }));
        } else if (type === 'energyDensityGamma') {
            // energyDensityGamma：与 energyDensity 同源，仅在 y 上做 gamma=0.1
          const gammaCorrect = (y) => {
            if (typeof y !== 'number' || !Number.isFinite(y)) return 0;
            if (y === 0) return 0;
            const s = y < 0 ? -1 : 1;
            const v = Math.pow(Math.abs(y), 0.1);
            return Number.isFinite(v) ? (s * v) : 0;
          };
            const baseKey = (slotConfig.orbital || '').replace(/[xyz]/g, '').replace(/_.*/, '');
            const energies = window.SlaterBasis && window.SlaterBasis[atomType]
              ? window.SlaterBasis[atomType].energies
              : null;
            const epsilon = (energies && energies[baseKey] !== undefined) ? energies[baseKey] : -0.5;
            theoryData = centers.map(r => {
              const y = epsilon * window.Hydrogen.radialPDF(orbitalParams.n, orbitalParams.l, r, atomZ, 1, atomType);
              return { x: r, y: gammaCorrect(y) };
            });
        } else if (type === 'energyCumulative') {
            // 能量累计理论曲线：与采样侧严格同定义：E(R)=∫ ε·P(r)dr
            const baseKey = (slotConfig.orbital || '').replace(/[xyz]/g, '').replace(/_.*/, '');
            const energies = window.SlaterBasis && window.SlaterBasis[atomType]
              ? window.SlaterBasis[atomType].energies
              : null;
            const epsilon = (energies && energies[baseKey] !== undefined) ? energies[baseKey] : -0.5;
            const dr = (hist && hist.edges && hist.edges.length >= 2)
              ? (hist.edges[1] - hist.edges[0])
              : (centers.length >= 2 ? (centers[1] - centers[0]) : 0.01);
            let cumSum = 0;
            theoryData = centers.map(r => {
              const density = epsilon * window.Hydrogen.radialPDF(orbitalParams.n, orbitalParams.l, r, atomZ, 1, atomType);
              cumSum += density * dr;
              return { x: r, y: cumSum };
            });
        } else if (type === 'energyCumulativeGamma') {
            // energyCumulativeGamma：理论 E 与 energyCumulative 同源，仅在 y 上做 gamma=0.1
          const gammaCorrect = (y) => {
            if (typeof y !== 'number' || !Number.isFinite(y)) return 0;
            if (y === 0) return 0;
            const s = y < 0 ? -1 : 1;
            const v = Math.pow(Math.abs(y), 0.1);
            return Number.isFinite(v) ? (s * v) : 0;
          };
            const baseKey = (slotConfig.orbital || '').replace(/[xyz]/g, '').replace(/_.*/, '');
            const energies = window.SlaterBasis && window.SlaterBasis[atomType]
              ? window.SlaterBasis[atomType].energies
              : null;
            const epsilon = (energies && energies[baseKey] !== undefined) ? energies[baseKey] : -0.5;
            const dr = (hist && hist.edges && hist.edges.length >= 2)
              ? (hist.edges[1] - hist.edges[0])
              : (centers.length >= 2 ? (centers[1] - centers[0]) : 0.01);
            let cumSum = 0;
            theoryData = centers.map(r => {
              const density = epsilon * window.Hydrogen.radialPDF(orbitalParams.n, orbitalParams.l, r, atomZ, 1, atomType);
              cumSum += density * dr;
              return { x: r, y: gammaCorrect(cumSum) };
            });
        }
        else if (type === 'localEnergy') {
          // compare-localEnergy：全部为理论派生量（无采样曲线）
          const baseKey = (slotConfig.orbital || '').replace(/[xyz]/g, '').replace(/_.*/, '');
          const energies = window.SlaterBasis && window.SlaterBasis[atomType]
            ? window.SlaterBasis[atomType].energies
            : null;
          const epsilon = (energies && energies[baseKey] !== undefined) ? energies[baseKey] : -0.5;

          const epsDensity = centers.map(r => {
            const P = window.Hydrogen.radialPDF(orbitalParams.n, orbitalParams.l, r, atomZ, 1, atomType);
            return { x: r, y: epsilon * P };
          });

          const Tdensity = centers.map(r => {
            const P = window.Hydrogen.radialPDF(orbitalParams.n, orbitalParams.l, r, atomZ, 1, atomType);
            const Tloc = (window.Hydrogen.calculateLocalKineticEnergy
              ? window.Hydrogen.calculateLocalKineticEnergy(r, epsilon, atomZ, atomType, baseKey)
              : window.Hydrogen.calculateReconstructedKineticEnergy(r, epsilon, atomZ, atomType, baseKey));
            return { x: r, y: Tloc * P };
          });

          datasets.push({
            label: `${displayName} ε·P(r)`,
            data: epsDensity,
            borderColor: `rgba(${colorValues.join(',')}, 0.95)`,
            backgroundColor: 'transparent',
            borderWidth: 2,
            pointRadius: 0,
            pointHoverRadius: 0,
            fill: false,
            tension: 0.2,
            borderDash: [6, 3],
            hidden: false
          });

          datasets.push({
            label: `${displayName} (ε - V_eff)·P(r)`,
            data: Tdensity,
            borderColor: `rgba(${colorValues.join(',')}, 0.85)`,
            backgroundColor: 'transparent',
            borderWidth: 2,
            pointRadius: 0,
            pointHoverRadius: 0,
            fill: false,
            tension: 0.2,
            borderDash: [2, 2],
            hidden: false
          });

          // localEnergy 已在此处直接 push 了两条曲线，不走后续统一的“理论曲线”逻辑
          theoryData = null;
        }
        else {
          // φ 方位角理论曲线
          theoryData = centers.map(phi => ({
            x: phi,
            y: window.Hydrogen.angularPDF_Phi(orbitalParams.angKey.m, orbitalParams.angKey.t, phi)
          }));
        }

        if (theoryData) {
          datasets.push({
            label: `${displayName} 理论`,
            data: theoryData,
            borderColor: `rgba(${colorValues.join(',')}, 0.9)`,
            backgroundColor: 'transparent',
            borderWidth: 2,
            pointRadius: 0,
            pointHoverRadius: 0,
            fill: false,
            tension: 0.2,
            borderDash: [6, 3], // 虚线样式
            // compare 的能量图只有理论曲线，必须默认可见，否则看起来像“被毁了”
            hidden: (type === 'potential' || type === 'dEdr') ? false : state.compareTheoryHidden
          });
        }
      }
    }

    applyThinLineStyle(datasets);
    state.chart.data.datasets = datasets;
    applyRollingAnimationPolicy(state.chart);
    applyRollingYAxisLock(state.chart, datasets, `compare:${type}`, type);
    console.log('开始更新图表，轨道数量:', Object.keys(orbitalDataMap).length);
    try {
      state.chart.update('none');
      console.log('图表更新完成');
    } catch (error) {
      console.error('图表更新失败:', error);
    }
  }

  function reset() {
    if (state.chart) {
      state.chart.data.labels = [];
      state.chart.data.datasets = [];
      state.chart.update();
    }
  }

  // 强制销毁并重建：用于“切换模式”时避免旧图表类型/样式残留
  function destroyChart() {
    if (state.chart) {
      try { state.chart.destroy(); } catch (e) { console.warn('图表销毁失败:', e); }
      state.chart = null;
    }
  }

  // 供外部调用的 API
  window.DataPanel = {
    init,
    reset,
    destroyChart,
    renderChartRadial,
    renderChartAngular,
    renderChartPhi,
    renderChartEnergyDensity,
    renderChartEnergyDensityGamma,
    renderChartEnergyCumulative,
    renderChartEnergyCumulativeGamma,
    renderChartEnergy, // 兼容
    renderChartPotential,
    renderChartDEdr,
    renderChartDEdrLog,
    renderChartLocalEnergy,
    renderChartCompare,
    state,
  };

  // 自动初始化
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
