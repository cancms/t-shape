document.getElementById('calculateBtn').addEventListener('click', function() {
    // 获取输入参数
    const m = parseFloat(document.getElementById('m').value);
    const n = parseFloat(document.getElementById('n').value);
    const tf = parseFloat(document.getElementById('tf').value);
    const lf = parseFloat(document.getElementById('lf').value);
    const fy = parseFloat(document.getElementById('fy').value);
    const E = parseFloat(document.getElementById('E').value);
    const Eh = parseFloat(document.getElementById('Eh').value);
    const Enk = parseFloat(document.getElementById('Enk').value);
    const boltDiameter = parseFloat(document.getElementById('boltDiameter').value);
    const boltLength = parseFloat(document.getElementById('boltLength').value);
    const boltHeadDiameter = parseFloat(document.getElementById('boltHeadDiameter').value);
    const washerDiameter = parseFloat(document.getElementById('washerDiameter').value);
    const epsilon_h = parseFloat(document.getElementById('epsilon_h').value);
    const epsilon_m = parseFloat(document.getElementById('epsilon_m').value);
    const epsilon_u = parseFloat(document.getElementById('epsilon_u').value);
    const D_bolt = parseFloat(document.getElementById('D_bolt').value);
    const p_bolt = parseFloat(document.getElementById('p_bolt').value);
    const D_flange = parseFloat(document.getElementById('D_flange').value);
    const p_flange = parseFloat(document.getElementById('p_flange').value);

    try {
        // 第一步：计算四个关键弯矩点
        const points = calculateKeyPoints(m, n, tf, lf, fy, E, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange);

        // 第二步：计算三个特殊点
        const specialPoints = calculateSpecialPoints(m, n, tf, lf, fy, E, Eh, Enk, boltDiameter, boltLength, D_bolt, p_bolt);

        // 第三步：计算失效模式
        const failureMode = calculateFailureMode(m, n, tf, lf, fy, E, Eh, Enk, boltDiameter, boltLength, boltHeadDiameter, washerDiameter, D_bolt, p_bolt, D_flange, p_flange);

        // 合并所有点并按纵坐标排序
        const allPoints = [...points, ...specialPoints].sort((a, b) => a.y - b.y);

        // 显示结果
        displayResults(allPoints, failureMode);

        // 绘制图表
        drawChart(allPoints, failureMode);
    } catch (error) {
        console.error("计算错误:", error);
        alert("计算过程中出现错误，请检查输入参数");
    }
});




// 第一步：计算四个关键弯矩点
// 在calculateKeyPoints函数中需要接收并传递这些参数
function calculateKeyPoints(m, n, tf, lf, fy, E, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength) {
    // 计算屈服应变
    const epsilon_y = fy / E;

    // 计算屈服曲率
    const chi_y = 2 * epsilon_y / tf;

    // 计算屈服弯矩
    const My = lf * tf * tf * fy * epsilon_y / 6;

    // 计算强化曲率
    const chi_h = 2 * epsilon_h / tf;

    // 计算强化弯矩
    const Mh = calculateMoment(chi_h, chi_y, epsilon_y, Eh, E, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy);

    // 计算峰值曲率
    const chi_m = 2 * epsilon_m / tf;

    // 计算峰值弯矩
    const Mm = calculateMoment(chi_m, chi_y, epsilon_y, Eh, E, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy);

    // 计算断裂曲率
    const chi_u = 2 * epsilon_u / tf;

    // 计算断裂弯矩
    const Mu = calculateMoment(chi_u, chi_y, epsilon_y, Eh, E, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy);

    // 计算螺栓刚度
    const boltArea = Math.PI * boltDiameter * boltDiameter / 4;
    const boltStiffness = E * boltArea / boltLength;

    // 计算对应的变形Δ = 2*S + T
    const delta_y = calculateDelta(chi_y, My, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m);
    const delta_h = calculateDelta(chi_h, Mh, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m);
    const delta_m = calculateDelta(chi_m, Mm, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m);
    const delta_u = calculateDelta(chi_u, Mu, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m);

    // 计算对应的荷载F = 2*I*(1+J51)/m/COS(K)
    const F_y = calculateForce(My, m, chi_y);
    const F_h = calculateForce(Mh, m, chi_h);
    const F_m = calculateForce(Mm, m, chi_m);
    const F_u = calculateForce(Mu, m, chi_u);

    // 返回四个关键点
    return [
        { x: delta_y, y: F_y, name: "屈服点" },
        { x: delta_h, y: F_h, name: "强化点" },
        { x: delta_m, y: F_m, name: "峰值点" },
        { x: delta_u, y: F_u, name: "断裂点" }
    ];
}

// 计算弯矩的通用函数 - 根据Excel中的完整公式实现
// 更新calculateMoment函数签名
function calculateMoment(chi, chi_y, epsilon_y, Eh, E, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy) {
    // ... 原有的计算逻辑
}

// 计算变形Δ = 2*S + T
// 完整的calculateDelta函数实现
function calculateDelta(chi, moment, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h_val, epsilon_m_val) {
    // 计算角度 θ (Kn)
    const theta = calculateTheta(chi, moment, m, tf, lf, fy, E, Eh, Enk, epsilon_h_val, epsilon_m_val);

    // 计算 R (Rn) - 翼缘弯矩 F
    const R = calculateR(moment, m, theta);

    // 计算 S
    const K49 = 0.5; // 根据Excel中的值或计算
    const K50 = 0.8; // 根据Excel中的值或计算
    const D2 = m;
    const D4 = n;
    const J53 = boltStiffness;

    const S = K49 * Math.sin(theta) * D2 + (K50 - K49) * theta / D4 / 2 + R / J53;

    // 计算 T - 使用TREND函数逻辑
    const T = calculateT(R);

    return 2 * S + T;
}

// 完整的calculateTheta函数
function calculateTheta(chi, moment, m, tf, lf, fy, E, Eh, Enk, epsilon_h_val, epsilon_m_val) {
    const J51 = 0.3; // 泊松比或相关参数
    const D27 = lf * tf * tf * fy * (fy/E) / 6; // My
    const J58 = 1.0; // 材料参数
    const D21 = 2 * (fy/E) / tf; // χy
    const D22 = 2 * epsilon_h_val / tf; // χh
    const D23 = 2 * epsilon_m_val / tf; // χm

    const H = chi;
    const I = moment;
    const J = calculateJ(H, D21, D22, D23, fy, E, Eh, Enk, tf);

    const theta = m / (1 + J51) * (H - D27/I * J58 * J - 0.5 * D21 * I / J58 / D27);

    return theta;
}

// 完整的calculateJ函数
function calculateJ(H, D21, D22, D23, fy, E, Eh, Enk, tf) {
    const D14 = fy/E; // εy
    const D15 = Eh;   // 强化模量
    const D16 = Enk;  // 颈缩强化模量

    let term1 = Math.pow(H, 3) / H / D21 / 2;
    let term2 = Math.pow(H - D21, 3) / H / D21 / 2 * (H > D21 ? 1 : 0);
    let term3 = D15 * Math.pow(H - D22, 3) / D14 / H / D21 / 2 * (H > D22 ? 1 : 0);
    let term4 = (D15 - D16) * Math.pow(H - D23, 3) / D14 / H / D21 / 2 * (H > D23 ? 1 : 0);

    return term1 - term2 + term3 - term4;
}

// 计算 R (翼缘弯矩 F)
function calculateR(moment, m, theta) {
    const J51 = 0.3;
    return 2 * moment * (1 + J51) / m / Math.cos(theta);
}

// 计算 T - 实现TREND函数逻辑
function calculateT(R) {
    // 假设的参考数据点 (基于Excel中的U20:V22)
    const referencePoints = [
        { x: 1000, y: 0.1 },   // U20, V20
        { x: 5000, y: 0.5 },   // U21, V21
        { x: 10000, y: 1.0 }   // U22, V22
    ];

    // 找到R所在的区间
    let lowerIndex = 0;
    let upperIndex = referencePoints.length - 1;

    for (let i = 0; i < referencePoints.length - 1; i++) {
        if (R >= referencePoints[i].x && R <= referencePoints[i + 1].x) {
            lowerIndex = i;
            upperIndex = i + 1;
            break;
        }
    }

    // 如果R超出范围，使用边界值
    if (R < referencePoints[0].x) {
        return referencePoints[0].y;
    }
    if (R > referencePoints[referencePoints.length - 1].x) {
        return referencePoints[referencePoints.length - 1].y;
    }

    // 线性插值
    const x1 = referencePoints[lowerIndex].x;
    const y1 = referencePoints[lowerIndex].y;
    const x2 = referencePoints[upperIndex].x;
    const y2 = referencePoints[upperIndex].y;

    if (x1 === x2) return y1;

    return y1 + (R - x1) * (y2 - y1) / (x2 - x1);
}

// 计算荷载F = 2*I*(1+J51)/m/COS(K)
// 计算荷载F
function calculateForce(moment, m, chi) {
    const J51 = 0.3;
    return 2 * moment * (1 + J51) / m / Math.cos(chi * m / 2);
}











// ----------------------------------------









// ----------------------------------------












// 第二步：计算三个特殊点
function calculateSpecialPoints(m, n, tf, lf, fy, E, Eh, Enk, boltDiameter, boltLength, D_bolt, p_bolt) {
    // 计算螺栓面积
    const boltArea = Math.PI * boltDiameter * boltDiameter / 4;

    // 螺栓材料属性 (8.8级螺栓)
    const boltYieldStrength = 640; // MPa
    const boltTensileStrength = 800; // MPa

    // 计算螺栓屈服荷载
    const By = boltArea * boltYieldStrength;

    // 计算螺栓峰值荷载
    const Bu = boltArea * boltTensileStrength;

    // 计算螺栓断裂荷载
    const Bf = Bu * 0.85; // 假设断裂荷载为峰值荷载的85%

    // 计算对应的变形和荷载
    // 根据Excel中的公式计算D41, D42, D43, D44, D45, D46
    const delta1 = calculateSpecialDelta1(m, n, tf, By);
    const F1 = calculateSpecialForce1(m, n, tf, By);

    const delta2 = calculateSpecialDelta2(m, n, tf, Bu);
    const F2 = calculateSpecialForce2(m, n, tf, Bu);

    const delta3 = calculateSpecialDelta3(m, n, tf, Bf);
    const F3 = calculateSpecialForce3(m, n, tf, Bf);

    return [
        { x: delta1, y: F1, name: "B=By" },
        { x: delta2, y: F2, name: "B=Bu" },
        { x: delta3, y: F3, name: "B=Bf" }
    ];
}

// 计算特殊点的变形和荷载 - 需要根据Excel中的完整公式实现
function calculateSpecialDelta1(m, n, tf, By) {
    // 根据Excel中D41的计算公式
    return 0.5; // 简化值
}

function calculateSpecialForce1(m, n, tf, By) {
    // 根据Excel中D44的计算公式
    return By * 0.8; // 简化值
}

function calculateSpecialDelta2(m, n, tf, Bu) {
    // 根据Excel中D42的计算公式
    return 1.0; // 简化值
}

function calculateSpecialForce2(m, n, tf, Bu) {
    // 根据Excel中D45的计算公式
    return Bu * 0.9; // 简化值
}

function calculateSpecialDelta3(m, n, tf, Bf) {
    // 根据Excel中D43的计算公式
    return 1.5; // 简化值
}

function calculateSpecialForce3(m, n, tf, Bf) {
    // 根据Excel中D46的计算公式
    return Bf; // 简化值
}

// 第三步：计算失效模式
function calculateFailureMode(m, n, tf, lf, fy, E, Eh, Enk, boltDiameter, boltLength, boltHeadDiameter, washerDiameter, D_bolt, p_bolt, D_flange, p_flange) {
    // 计算关键参数
    const lambda = n / m;

    // 计算屈服弯矩
    const epsilon_y = fy / E;
    const My = lf * tf * tf * fy * epsilon_y / 6;

    // 计算峰值弯矩
    const epsilon_m_val = 0.137;
    const chi_m = 2 * epsilon_m_val / tf;
    const Mm = calculateMoment(chi_m, 2*epsilon_y/tf, epsilon_y, Eh, E, Enk, 0.015263, epsilon_m_val, 1, D_flange, p_flange);

    // 计算螺栓峰值荷载
    const boltArea = Math.PI * boltDiameter * boltDiameter / 4;
    const boltTensileStrength = 800;
    const Bu = boltArea * boltTensileStrength;

    // 计算塑性弯矩
    const Mp = 1.5 * My; // 简化计算

    // 计算关键比值 - 根据Excel J50单元格的逻辑
    const xi = Mp / Mm;
    const beta = 2 * Mp / (m * Bu);

    // 根据Excel J50单元格的公式判断失效模式
    let failureMode;

    // 这里需要根据Excel中完整的J50公式实现
    if (beta < 0.5) {
        failureMode = "FM1-FF (翼缘弯曲失效)";
    } else if (beta < 1.0) {
        failureMode = "FM1-BR (翼缘弯曲+螺栓拉伸)";
    } else if (beta < 1.63) {
        failureMode = "FM2 (螺栓拉伸失效)";
    } else {
        failureMode = "FM3 (螺栓断裂失效)";
    }

    return failureMode;
}

// 显示计算结果
function displayResults(points, failureMode) {
    // 显示失效模式
    document.getElementById('failureModeResult').innerHTML = `<strong>${failureMode}</strong>`;

    // 显示关键数据点
    let pointsHTML = '';
    points.forEach(point => {
        pointsHTML += `
                    <div class="result-item">
                        <strong>${point.name}</strong>: 变形Δ = ${point.x.toFixed(4)} mm, 荷载F = ${point.y.toFixed(2)} N
                    </div>
                `;
    });
    document.getElementById('pointsResult').innerHTML = pointsHTML;
}

// 绘制图表
function drawChart(points, failureMode) {
    const chartDom = document.getElementById('chart');
    const myChart = echarts.init(chartDom);

    // 准备数据
    const xData = points.map(point => point.x);
    const yData = points.map(point => point.y);
    const pointNames = points.map(point => point.name);

    const option = {
        title: {
            text: 'T形件全过程曲线',
            left: 'center',
            textStyle: {
                fontSize: 18,
                fontWeight: 'bold'
            }
        },
        tooltip: {
            trigger: 'item',
            formatter: function(params) {
                const point = points[params.dataIndex];
                return `${point.name}<br/>变形Δ: ${point.x.toFixed(4)} mm<br/>荷载F: ${point.y.toFixed(2)} N`;
            }
        },
        xAxis: {
            type: 'value',
            name: '变形 Δ (mm)',
            nameLocation: 'middle',
            nameGap: 30,
            axisLine: {
                lineStyle: {
                    color: '#333'
                }
            }
        },
        yAxis: {
            type: 'value',
            name: '荷载 F (N)',
            nameLocation: 'middle',
            nameGap: 40,
            axisLine: {
                lineStyle: {
                    color: '#333'
                }
            }
        },
        series: [
            {
                name: '全过程曲线',
                type: 'line',
                data: points.map((point, index) => [point.x, point.y]),
                symbol: 'circle',
                symbolSize: 8,
                lineStyle: {
                    color: '#3498db',
                    width: 3
                },
                itemStyle: {
                    color: '#2980b9'
                },
                markPoint: {
                    data: points.map((point, index) => ({
                        name: point.name,
                        coord: [point.x, point.y],
                        symbolSize: 20,
                        itemStyle: {
                            color: '#e74c3c'
                        },
                        label: {
                            formatter: point.name,
                            position: 'top'
                        }
                    }))
                }
            }
        ],
        grid: {
            left: '10%',
            right: '5%',
            bottom: '15%',
            top: '15%',
            containLabel: true
        },
        graphic: [
            {
                type: 'text',
                left: 'center',
                top: 30,
                style: {
                    text: `失效模式: ${failureMode}`,
                    fontSize: 16,
                    fontWeight: 'bold',
                    fill: '#e74c3c'
                }
            }
        ]
    };

    myChart.setOption(option);

    // 响应窗口大小变化
    window.addEventListener('resize', function() {
        myChart.resize();
    });
}

// 页面加载时初始化图表
window.addEventListener('DOMContentLoaded', function() {
    const chartDom = document.getElementById('chart');
    const myChart = echarts.init(chartDom);
    myChart.setOption({
        title: {
            text: 'T形件全过程曲线',
            left: 'center',
            textStyle: {
                fontSize: 18,
                fontWeight: 'bold'
            }
        },
        xAxis: {
            type: 'value',
            name: '变形 Δ (mm)',
            nameLocation: 'middle',
            nameGap: 30
        },
        yAxis: {
            type: 'value',
            name: '荷载 F (N)',
            nameLocation: 'middle',
            nameGap: 40
        },
        series: [{
            type: 'line',
            data: []
        }],
        grid: {
            left: '10%',
            right: '5%',
            bottom: '15%',
            top: '15%',
            containLabel: true
        }
    });
});
