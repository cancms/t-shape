// 全部关键点弯矩 ： 屈服弯矩, 强化弯矩，峰值弯矩，断裂弯矩， 塑性弯矩
let allMoments = {
    My: 0, // 屈服弯矩
    Mh: 0, // 强化弯矩
    Mm: 0, // 峰值弯矩
    Mu: 0, // 断裂弯矩
    Mp: 0, // 塑性弯矩
};


// 位移相容系数k49 , k50
let k49k50Params = {};

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
        // const points = calculateKeyPoints(m, n, tf, lf, fy, E, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange);
        const points = calculateKeyPoints(m, n, tf, lf, fy, E, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);

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
    // const My = lf * tf * tf * fy * epsilon_y / 6;

    // const My = lf * tf * tf * fy  / 6;

    // console.log("屈服弯矩0", My)

    // 计算强化曲率
    const chi_h = 2 * epsilon_h / tf;


    allMoments = calculateAllMoments(m, n, tf, lf, fy, E, Eh, Enk, epsilon_h, epsilon_m, epsilon_u);

    console.log("allMoments", allMoments)

    // 屈服弯矩
    const My = allMoments.My;


    // 计算强化弯矩
    // const Mh = calculateMoment(chi_h, chi_y, epsilon_y, Eh, E, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy);
    const Mh = allMoments.Mh;

    // 计算峰值曲率
    const chi_m = 2 * epsilon_m / tf;

    // 计算峰值弯矩
    // const Mm = calculateMoment(chi_m, chi_y, epsilon_y, Eh, E, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy);
    const Mm = allMoments.Mm;



    // 计算断裂曲率
    const chi_u = 2 * epsilon_u / tf;

    // 计算断裂弯矩
    // const Mu = calculateMoment(chi_u, chi_y, epsilon_y, Eh, E, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy);
    const Mu = allMoments.Mu;

    // 计算螺栓刚度
    // const boltArea = Math.PI * boltDiameter * boltDiameter / 4;
    // const boltStiffness = E * boltArea / boltLength;
    const boltStiffness = 0; // 用单元格J53的值

    // 计算对应的变形Δ = 2*S + T
    // const delta_y = calculateDelta(chi_y, My, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);
    // const delta_h = calculateDelta(chi_h, Mh, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);
    // const delta_m = calculateDelta(chi_m, Mm, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);
    // const delta_u = calculateDelta(chi_u, Mu, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);

    const the4KeyPoints = [];
    let thePoint = {};




    // My, 屈服弯矩点， 屈服点，   / 1000 ， 用千牛作为纵坐标的单位
    thePoint = calculateDelta(chi_y, My, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);
    the4KeyPoints.push({ x: thePoint.U, y: thePoint.R  / 1000, name: "屈服点" });

    // Mh, 屈服弯矩点， 强化点
    thePoint = calculateDelta(chi_h, Mh, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);
    the4KeyPoints.push({ x: thePoint.U, y: thePoint.R  / 1000, name: "强化点" });

    // Mm, 峰值弯矩点， 峰值点
    thePoint = calculateDelta(chi_m, Mm, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);
    the4KeyPoints.push({ x: thePoint.U, y: thePoint.R  / 1000, name: "峰值点" });

    // Mu, 断裂弯矩点， 断裂点
    thePoint = calculateDelta(chi_u, Mu, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h, epsilon_m, epsilon_u, D_flange, p_flange, boltDiameter, boltLength);
    the4KeyPoints.push({ x: thePoint.U, y: thePoint.R  / 1000, name: "断裂点" });

    // 计算对应的荷载F = 2*I*(1+J51)/m/COS(K)
    // const F_y = calculateForce(My, m, chi_y);
    // const F_h = calculateForce(Mh, m, chi_h);
    // const F_m = calculateForce(Mm, m, chi_m);
    // const F_u = calculateForce(Mu, m, chi_u);

    // 返回四个关键点
    return the4KeyPoints;
    // return [
    //     { x: delta_y, y: F_y, name: "屈服点" },
    //     { x: delta_h, y: F_h, name: "强化点" },
    //     { x: delta_m, y: F_m, name: "峰值点" },
    //     { x: delta_u, y: F_u, name: "断裂点" }
    // ];
}


// 屈服弯矩, lf 宽度， tf厚度， fy(屈服应力410) = E(弹性模量) * epsilon_y(屈服应变)
function calculateMy(tf, lf, fy) {
    const My = lf * tf * tf * fy  / 6;
    return My; // 屈服弯矩
}

function calculateAllMoments(m, n, tf, lf, fy, E, Eh, Enk,
                             epsilon_h = 0.015263, epsilon_m = 0.137, epsilon_u = 1) {

    // 强化弯矩
    function calculateMh(chi_h, chi_y, My) {
        // chi_h = 强化曲率 = 2 * εh / tf
        // chi_y = 屈服曲率 = 2 * εy / tf
        // My = 屈服弯矩

        const ratio = chi_h / chi_y;
        const term1 = ratio;
        const term2 = 0.5 * (3 - 2 * ratio - Math.pow(chi_y / chi_h, 2));

        return (term1 + term2) * My;
    }

    // 峰值弯矩
    function calculateMm(chi_m, chi_y, chi_h, Eh, E, epsilon_y, My) {
        // chi_m = 峰值曲率 = 2 * εm / tf
        // chi_h = 强化曲率 = 2 * εh / tf
        // Eh = 强化模量
        // E = 弹性模量
        // epsilon_y = 屈服应变 = fy/E

        const ratio_m = chi_m / chi_y;
        const ratio_h = chi_h / chi_y;

        // 基础项
        const term1 = ratio_m;
        const term2 = 0.5 * (3 - 2 * ratio_m - Math.pow(chi_y / chi_m, 2));

        // 强化修正项
        const term3 = 0.5 * (Eh / E) * (ratio_m - ratio_h);
        const term4 = (1 - chi_h / chi_m) * (2 + chi_h / chi_m);

        return (term1 + term2 + term3 * term4) * My;
    }


    // 断裂弯矩
    function calculateMu(chi_u, chi_y, chi_h, chi_m, Eh, Enk, E, epsilon_y, My) {
        // chi_u = 断裂曲率 = 2 * εu / tf
        // chi_m = 峰值曲率 = 2 * εm / tf
        // Enk = 颈缩强化模量

        const ratio_u = chi_u / chi_y;
        const ratio_h = chi_h / chi_y;
        const ratio_m = chi_m / chi_y;

        // 基础项
        const term1 = ratio_u;
        const term2 = 0.5 * (3 - 2 * ratio_u - Math.pow(chi_y / chi_u, 2));

        // 强化修正项
        const term3 = 0.5 * (Eh / E) * (ratio_u - ratio_h);
        const term4 = (1 - chi_h / chi_u) * (2 + chi_h / chi_u);

        // 颈缩修正项
        const term5 = 0.5 * ((Eh - Enk) / E) * (ratio_u - ratio_m);
        const term6 = (1 - chi_m / chi_u) * (2 + chi_m / chi_u);

        return (term1 + term2 + term3 * term4 - term5 * term6) * My;
    }

    // 基本参数计算
    const epsilon_y = fy / E;
    const chi_y = 2 * epsilon_y / tf;
    const chi_h = 2 * epsilon_h / tf;
    const chi_m = 2 * epsilon_m / tf;
    const chi_u = 2 * epsilon_u / tf;

    // 屈服弯矩
    // const My = lf * tf * tf * fy * epsilon_y / 6;
    const My = calculateMy(tf, lf, fy);

    // 强化弯矩
    const Mh = calculateMh(chi_h, chi_y, My);

    // 峰值弯矩
    const Mm = calculateMm(chi_m, chi_y, chi_h, Eh, E, epsilon_y, My);

    // 断裂弯矩
    const Mu = calculateMu(chi_u, chi_y, chi_h, chi_m, Eh, Enk, E, epsilon_y, My);

    // 塑性弯矩 (Excel D31单元格)
    const Mp = My * 3 / 2;

    return {
        My: My,
        Mh: Mh,
        Mm: Mm,
        Mu: Mu,
        Mp: Mp,
        curvatures: {
            chi_y: chi_y,
            chi_h: chi_h,
            chi_m: chi_m,
            chi_u: chi_u
        }
    };
}




// 位移相容系数k49 , k50
// 从基础输入参数开始计算K49和K50
function calculateK49K50(m, n, tf, lf, fy, E, Eh, Enk, epsilon_h_val, epsilon_m_val, epsilon_u, D_flange, p_flange,  boltDiameter, boltLength) {
    // 从基础输入参数开始计算K49和K50
    function calculateK49K50FromInput(m, n, tf, lf, fy, E, Eh, Enk, boltDiameter, boltLength) {

        const D2 = m;

        // 第一步：计算基本几何参数
        const D4 = n / m;  // λ = n/m

        const D5 = tf;
        const D6 = lf;
        const D14 = E;

        const D18 = p_flange; // 翼缘率强化指数 (p)
        const D17 = D_flange;  // 翼缘率强化指数 (D)

        const J39 = 0;
        const J43 = 2 * 0.1 * D2

        // 第二步：计算关键弯矩值
        const epsilon_y = fy / E;

        // console.log("epsilon_y", epsilon_y)
        // const My = lf * tf * tf * fy * epsilon_y / 6;  // 屈服弯矩
        // const My = lf * tf * tf * fy  / 6;  // 屈服弯矩
        const My = allMoments.My;

        // console.log("屈服弯矩", My)

        // 计算峰值弯矩 (简化)
        // const epsilon_m = epsilon_m_val; // 0.137;
        // const chi_m = 2 * epsilon_m / tf;
        // const Mm = calculateMoment(chi_m, 2*epsilon_y/tf, epsilon_y, Eh, E, Enk, epsilon_h_val, epsilon_m, epsilon_u, D_flange, p_flange, tf, lf, fy);
        const Mm = allMoments.Mm;
        // console.log("计算峰值弯矩", Mm)


        // 断裂弯矩
        const Mu = allMoments.Mu;


        // 第三步：计算螺栓相关参数
        const boltArea = Math.PI * boltDiameter * boltDiameter / 4;
        const boltTensileStrength = 930; // C63
        const Bu = boltArea * boltTensileStrength; // 螺栓峰值荷载, BU(D45) = =D35^2*PI()/4*C63

        // 第四步：计算J系列参数
        // const Mp = 1.5 * My;  // 塑性弯矩 (简化)
        const Mp =  allMoments.Mp;  // 塑性弯矩

        // J40 - 需要根据Excel确定，这里假设一个值, =1+(J39/D36/D38)^(1/D39)
        const J40 = 1.0;



        // J46 --- ξ;
        // J47 --- β;
        // J48 --- β1;
        // J49 --- β2;

        // const J46 = Mp / Mm; // 这里用Mm代替Mu
        // J46: ξ = Mp/Mu，
        const J46 = Mp / Mu;

        // J42 = DM (弯矩放大系数) - 需要计算
        const J42 = calculateDM(m, tf, D_flange, p_flange);

        // J47 = β = 2Mp/(m*Bu)
        const J47 = 2 * Mp / (m * Bu);

        // J48, J49 - 失效模式阈值
        // J48 --- β1; J49 --- β2;
        const J48 = calculateJ48(J46, J42, D4);
        const J49 = calculateJ49Threshold(J40, J46, J42, D4);

        // 第五步：计算失效模式J50
        const J50 = calculateJ50(J47, J48, J49, J40);

        //  计算失效模式过渡系数
        const J51 = calculateJ51(J47, D4, J46, J42, J40);


        // 第八步：最终计算K49和K50
        const K49 = calculateK49(J47, J48, J49);
        // const K50 = calculateK50(J50, J47, J48, J49, N49, m, tf, D42);
        const K50 = calculateK50(J47, J40);

        // 螺栓刚度：(螺栓材料的强度)
        const J53 = calculateJ53(D14, D6, D5, D2);

        // 翼缘速率 mm/s
        const J54 = calculateJ54(J39, J47, J49, K49, D4);

        // 翼缘边缘应变率
        const J57 = calculateJ57(J54, D5, J43, D2);

        // 率强化修正系数
        const J58 = calculateJ58(D18, J57, D17);





        return {
            K49: K49,
            K50: K50,
            intermediateValues: {
                D4: D4,
                J42: J42,
                J46: J46,
                J47: J47,
                J48: J48,
                J49: J49,
                J50: J50,
                J51: J51,
                J53: J53,
                J58:J58,


            }
        };
    }

// 辅助函数实现

    function calculateDM(m, tf, D_flange, p_flange) {

        const D2 = m;
        const D5 = tf;
        const J43 = 2 * 0.1 * D2;
        const J39 = 0;

        const J41 = J39/2 * D5 / J43 / D2;
        const D17 = D_flange;

        // 弯矩放大系数计算
        // let strainRate = 0.001;
        const strainRate = J41 / D17;
        // const D18 = 0.1; // 率强化指数 (假设)
        const D18 = p_flange; // 翼缘率强化指数 (p)
        return 1 + 2 * D18 / (1 + 2 * D18) * Math.pow(strainRate, 1/D18);
    }

    // J48 --- β1;
    function calculateJ48(J46, J42, D4) {
        // J48阈值计算
        return 2 * D4 * J46 / ((1 + 2 * D4) * J42);
    }


    // = J40*J46 / (J42*J40-J46)  *  (SQRT(1+4*D4*(J42*J40-J46)/J46/(1+2*D4)*J40)-1)
    function calculateJ49Threshold(J40, J46, J42, D4) {
        // J49阈值计算
        const numerator = J40 * J46;
        const denominator = J42 * J40 - J46;

        if (denominator === 0) return 2 * J40;

        const sqrtPart = 1 + 4 * D4 * (J42 * J40 - J46) / J46 / (1 + 2 * D4) * J40;
        const trendValue = numerator / denominator * (Math.sqrt(Math.max(sqrtPart, 0)) - 1);

        return Math.min(trendValue, 2 * J40);
    }

    // 计算失效模式J50
    // =IF(J47-J48>0,IF(J47-J49>0,IF(J47-2*J40>0,"FM3","FM2"),"FM1-BR"),"FM1-FF")
    function calculateJ50(J47, J48, J49, J40) {
        if (J47 - J48 > 0) {
            if (J47 - J49 > 0) {
                if (J47 - 2 * J40 > 0) {
                    return "FM3";
                } else {
                    return "FM2";
                }
            } else {
                return "FM1-BR";
            }
        } else {
            return "FM1-FF";
        }
    }


    // 物理意义说明
    // J51是一个失效模式过渡系数，其计算逻辑：
    // J47 < 阈值：返回1，表示完全弹性状态
    // J47 > 2×J40：返回0，表示完全塑性状态
    // 中间范围：使用复杂分式计算过渡系数，反映几何参数和材料参数的综合影响
    // 这个系数在变形计算中起到重要的修正作用，确保不同受力状态下的计算准确性。
    function calculateJ51(J47, D4, J46, J42, J40) {
        // Excel公式:
        // =IF(J47<2*D4*J46/(1+2*D4)/J42,1,IF(J47>2*J40,0,D4*(2*J40-J47)/(((1+2*D4)/J46*J40*J42-D4)*J47)))

        // 计算第一个条件阈值
        const threshold1 = (2 * D4 * J46) / ((1 + 2 * D4) * J42);

        if (J47 < threshold1) {
            return 1;
        } else if (J47 > 2 * J40) {
            return 0;
        } else {
            // 计算复杂的分式
            const numerator = D4 * (2 * J40 - J47);
            const denominatorTerm1 = ((1 + 2 * D4) / J46) * J40 * J42;
            const denominator = (denominatorTerm1 - D4) * J47;

            // 避免除零错误
            if (Math.abs(denominator) < 1e-10) {
                return 0.5; // 返回中间值作为默认
            }

            return numerator / denominator;
        }
    }


    // 根据Excel上下文，这个公式计算的是螺栓刚度：(螺栓材料的强度)
    // D14：屈服强度 fy (MPa)
    // D6：有效宽度 lf (mm)
    // D5：翼缘厚度 tf (mm)
    // D2：螺栓到翼缘边缘距离 m (mm)
    function calculateJ53(D14, D6, D5, D2) {
        // Excel公式: =0.85*D14*D6*D5^3/D2^3

        // 避免除零错误
        if (D2 === 0) {
            return 0; // 如果D2为0，返回0
        }

        // 直接计算
        return 0.85 * D14 * D6 * Math.pow(D5, 3) / Math.pow(D2, 3);
    }



    // 根据Excel上下文，这个公式计算的是转动速度： (翼缘速率 mm/s)
    // J39：加载速率 v (mm/s)
    // J47：β值 = 2Mp/(m×Bu)
    // J49：失效模式阈值
    // K49：失效模式过渡系数
    // D4：几何参数 λ = n/m
    function calculateJ54(J39, J47, J49, K49, D4) {
        // Excel公式: =J39/2*IF(J47>J49,K49/(1+D4),1-D4*K49/(1+D4))

        // 第一部分: J39/2
        const baseTerm = J39 / 2;

        // 条件判断部分
        let conditionTerm;
        if (J47 > J49) {
            // J47 > J49 的情况: K49/(1+D4)
            conditionTerm = K49 / (1 + D4);
        } else {
            // J47 <= J49 的情况: 1 - D4*K49/(1+D4)
            conditionTerm = 1 - (D4 * K49) / (1 + D4);
        }

        // 最终结果
        return baseTerm * conditionTerm;
    }



    // 根据Excel中的上下文，这个公式计算的是翼缘边缘应变率：
    // J54：转动速度 θ· (rad/s)
    // D5：翼缘厚度 tf (mm)
    // J43：塑性铰长度 2lp (mm)
    // D2：螺栓到翼缘边缘距离 m (mm)
    function calculateJ57(J54, D5, J43, D2) {
        // Excel公式: =J54*D5/J43/D2

        // 避免除零错误
        if (J43 === 0 || D2 === 0) {
            return 0; // 如果分母为0，返回0
        }

        // 直接计算
        return (J54 * D5) / (J43 * D2);
    }


    // J58是一个率强化修正系数，其物理意义：
    // D18：率强化指数参数
    // J57：当前应变率或相关速率参数
    // D17：参考应变率或基准速率参数
    function calculateJ58(D18, J57, D17) {
        // Excel公式: =1+2*D18/(1+2*D18)*(J57/D17)^(1/D18)

        // 避免除零错误
        if (D18 === 0) {
            return 1; // 如果D18为0，返回基准值1
        }

        if (D17 === 0) {
            return 1; // 如果D17为0，避免除零错误
        }

        // 计算第一部分系数
        const coefficient = (2 * D18) / (1 + 2 * D18);

        // 计算幂的底数
        const base = J57 / D17;

        // 避免负数的分数次幂
        if (base < 0 && 1/D18 % 2 === 0) {
            return 1; // 如果底数为负且指数分母为偶数，返回基准值
        }

        // 计算幂
        const powerTerm = Math.pow(Math.abs(base), 1 / D18);

        // 如果底数为负且指数分母为奇数，处理符号
        let signedPowerTerm = base < 0 ? -powerTerm : powerTerm;

        // 最终计算
        return 1 + coefficient * signedPowerTerm;
    }



    function calculateD42(m, tf, boltDiameter) {
        // D42变形参数计算
        return 0.1 * m + 0.5 * tf + 0.2 * boltDiameter;
    }

    // =IF(J47>J49,0,IF(J47<J48,1,(J49-J47)/(J49-J48)))
    function calculateK49(J47, J48, J49) {
        if (J47 > J49) {
            return 0;
        } else if (J47 < J48) {
            return 1;
        } else {
            return (J49 - J47) / (J49 - J48);
        }
    }


    // =IF(J47>2*J40,0,IF(J47>1.63,(2-J47)/(2-1.63),1))
    // function calculateK50(J50, J47, J48, J49, N49, D2, D5, D42) {
    //     if (J50 === "FM1-BR") {
    //         const exponentBase = 0.33 * Math.pow(D2, 2) / (2 * D5 * D42);
    //         const exponent = 6 * Math.pow(exponentBase, 1.5);
    //         const ratio = (J47 - J48) / (J49 - J48);
    //
    //         return N49 + (1 - N49) * Math.pow(ratio, exponent);
    //     } else {
    //         return "--";
    //     }
    // }
    // K50单元格的计算方法
    function calculateK50(J47, J40) {
        // Excel公式: =IF(J47>2*J40,0,IF(J47>1.63,(2-J47)/(2-1.63),1))
        // J47 = β值 = 2Mp/(m*Bu)
        // J40 = 阈值参数

        if (J47 > 2 * J40) {
            return 0;
        } else if (J47 > 1.63) {
            return (2 - J47) / (2 - 1.63);
        } else {
            return 1;
        }
    }

    return calculateK49K50FromInput(m, n, tf, lf, fy, E, Eh, Enk, boltDiameter, boltLength);

}





// 计算变形Δ = 2*S + T
// 完整的calculateDelta函数实现
function calculateDelta(chi, moment, m, n, tf, lf, fy, E, boltStiffness, Eh, Enk, epsilon_h_val, epsilon_m_val,epsilon_u, D_flange, p_flange,  boltDiameter, boltLength) {


    // 计算 S
    // const K49 = 0.5; // 根据Excel中的值或计算
    // const K50 = 0.8; // 根据Excel中的值或计算
    k49k50Params = calculateK49K50(m, n, tf, lf, fy, E, Eh, Enk, epsilon_h_val, epsilon_m_val, epsilon_u, D_flange, p_flange,  boltDiameter, boltLength);
    const {K49, K50} = k49k50Params;

    console.log('k49k50Params', k49k50Params)

    // 计算角度 θ (Kn)
    const theta = calculateTheta(chi, moment, m, tf, lf, fy, E, Eh, Enk, epsilon_h_val, epsilon_m_val);

    // 计算 R (Rn) - 翼缘弯矩 F
    const R = calculateR(moment, m, theta);




    // x轴： 2Δ1+Δ2
    // S: Δ1    ;     T:  Δ2


    // 计算 T - 使用TREND函数逻辑
    // R(列) 为受力 F
    const T = calculateT(R);

    const D2 = m;

    // λ: = n/m
    const D4 = n/m;

    // 螺栓刚度：(螺栓材料的强度)
    // const J53 = boltStiffness;
    const J53 = k49k50Params.intermediateValues.J53;


    // S5=$K$49*SIN(K5)*$D$2+($K$50-$K$49)*T5/$D$4/2+R5/$J$53
    // const S = K49 * Math.sin(theta) * D2 + (K50 - K49) * theta / D4 / 2 + R / J53;
    const S = K49 * Math.sin(theta) * D2 + (K50 - K49) * T / D4 / 2 + R / J53;

    const U = 2 * S + T;

    // console.log("moment, theta, R(F) ", moment,  theta, R )
    console.log(`
        moment: ${moment},  
        theta: ${theta},  
        
        K49: ${K49},  
        K50: ${K50},  
        
        D2: ${D2},  
        D4: ${D4},  
        J53: ${J53},  
        
        R(F): ${R},  
        
        
        S: ${S},  
        T: ${T},  
        2Δ1+Δ2: ${U},  
    `);

    return {U, R};
    // return U;
}

// 完整的calculateTheta函数
function calculateTheta(chi, moment, m, tf, lf, fy, E, Eh, Enk, epsilon_h_val, epsilon_m_val) {
    // const J51 = 0.3; // 泊松比或相关参数
    const J51 = k49k50Params.intermediateValues.J51;

    // const D27 = lf * tf * tf * fy * (fy/E) / 6; // My
    // 屈服弯矩
    const D27 = allMoments.My

    const J58 = k49k50Params.intermediateValues.J58; // 率强化修正系数

    const D21 = 2 * (fy/E) / tf; // χy
    const D22 = 2 * epsilon_h_val / tf; // χh
    const D23 = 2 * epsilon_m_val / tf; // χm

    const H = chi;
    const I = moment;
    const J = calculateJ(H, D21, D22, D23, fy, E, Eh, Enk, tf);

    // console.log('moment, J column, H',moment, J, H)

    const theta = m / (1 + J51) * (H - D27/I * J58 * J - 0.5 * D21 * I / J58 / D27);

    return theta;
}

// 完整的calculateJ函数
function calculateJ(H, D21, D22, D23, fy, E, Eh, Enk, tf) {
    // const D14 = fy/E; // εy
    const D14 = E; // 弹性模量
    const D15 = Eh;   // 强化模量
    const D16 = Enk;  // 颈缩强化模量


    // console.log('H,  D21, D22, D23, D14, D15, D16', H, D21, D22, D23, E,  Eh, Enk);

    let term1 = Math.pow(H, 3) / H / D21 / 2;
    let term2 = Math.pow(H - D21, 3) / H / D21 / 2 * (H > D21 ? 1 : 0);
    let term3 = D15 * Math.pow(H - D22, 3) / D14 / H / D21 / 2 * (H > D22 ? 1 : 0);
    let term4 = (D15 - D16) * Math.pow(H - D23, 3) / D14 / H / D21 / 2 * (H > D23 ? 1 : 0);

    return term1 - term2 + term3 - term4;
}

// 计算 R (翼缘弯矩 F)
function calculateR(moment, m, theta) {
    const J51 = k49k50Params.intermediateValues.J51;
    return 2 * moment * (1 + J51) / m / Math.cos(theta);
}

// 计算 T - 实现TREND函数逻辑
// R(F): 为受力（R列的值）， 参考插值表的值，用trend（excel函数）插值算法根据Ma（弯矩的值）得到 T（Δ2）的值
function calculateT(R) {
    // 插值算法，示例1：简单线性回归（最小二乘法）
    // const tableYs = [100, 120, 150, 180, 200];
    // const tableYs = [1, 2, 3, 4, 5];
    // const inputXs = [6, 7, 8];
    // const predictions = trend(tableYs, tableXs, inputXs);
    // console.log('预测:', predictions);
    // 输出: [220, 240, 260] (近似值)
    function trend(y, x = [], new_x = [], constant = true) {
        // 如果x为空，使用默认序列 [1, 2, 3, ...]
        if (x.length === 0) {
            x = Array.from({length: y.length}, (_, i) => i + 1);
        }

        // 如果没有提供new_x，使用x作为预测点
        if (new_x.length === 0) {
            new_x = [...x];
        }

        const n = y.length;

        // 计算统计量
        const sumX = x.reduce((a, b) => a + b, 0);
        const sumY = y.reduce((a, b) => a + b, 0);
        const sumXY = x.reduce((sum, xi, i) => sum + xi * y[i], 0);
        const sumX2 = x.reduce((sum, xi) => sum + xi * xi, 0);

        let slope, intercept;

        if (constant) {
            // 有截距的线性回归: y = a + bx
            slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
            intercept = (sumY - slope * sumX) / n;
        } else {
            // 无截距的线性回归: y = bx
            slope = sumXY / sumX2;
            intercept = 0;
        }

        // 计算预测值
        return new_x.map(xi => intercept + slope * xi);
    }


    // 参考表, 假设的参考数据点 (基于Excel中的U21:V24)， 按 x 值从小到大排列
    // const referencePoints = [
    //     { x: 0, y: 0 },   // U21, V21
    //     { x: 32686, y: 0.07 },   // U22, V22
    //     { x: 37715, y: 3.95 },   // U24, V24
    //     { x: 46767, y: 1.30 },   // U23, V23
    // ];

    // excel表U20 - V22
    const referencePoints = [
        { x: 0, y: 0 },   // U20, V20
        { x: 0, y: 0 },   // U21, V21
        { x: 32686, y: 0.07 },   // U22, V22
        // { x: 37715, y: 3.95 },   // U24, V24
        // { x: 46767, y: 1.30 },   // U23, V23
    ];


    const tableXs = [];
    const tableYs = [];
    const inputXs = [R];
    referencePoints.map(p => {
        tableXs.push(p.x);
        tableYs.push(p.y);
    });
    const predictDelta2 = trend(tableYs, tableXs, inputXs);

    // T 列（Δ2）的值
    let T = 0;
    if (predictDelta2.length > 0 ) {
        T = predictDelta2[0];
    }
    return T;
}


// 计算荷载F - 改为千牛
function calculateForce(moment, m, chi) {
    const J51 = 0.3;
    // 除以1000将单位从牛(N)改为千牛(kN)
    return 2 * moment * (1 + J51) / m / Math.cos(chi * m / 2) / 1000;
}






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

// 在calculateSpecialForce函数中也做同样的修改
function calculateSpecialForce1(m, n, tf, By) {
    // 根据Excel中D44的计算公式，单位改为千牛
    return By * 0.8 / 1000; // 简化值，除以1000
}



function calculateSpecialDelta2(m, n, tf, Bu) {
    // 根据Excel中D42的计算公式
    return 1.0; // 简化值
}

function calculateSpecialForce2(m, n, tf, Bu) {
    // 根据Excel中D45的计算公式，单位改为千牛
    return Bu * 0.9 / 1000; // 简化值，除以1000
}

function calculateSpecialDelta3(m, n, tf, Bf) {
    // 根据Excel中D43的计算公式
    return 1.5; // 简化值
}

function calculateSpecialForce3(m, n, tf, Bf) {
    // 根据Excel中D46的计算公式，单位改为千牛
    return Bf / 1000; // 简化值，除以1000
}

// 第三步：计算失效模式
function calculateFailureMode() {
    const failureMode = k49k50Params.intermediateValues.J50 || 'FM1-FF';
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
                <strong>${point.name}</strong>: 变形Δ = ${point.x.toFixed(4)} mm, 荷载F = ${point.y.toFixed(2)} kN
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
                return `${point.name}<br/>变形Δ: ${point.x.toFixed(4)} mm<br/>荷载F: ${point.y.toFixed(2)} kN`;
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
            name: '荷载 F (kN)',  // 改为千牛
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
