

/**
 * @desc T形件全过程曲线计算器
 *    按m8p6的逻辑来算
 *
 *
 * @constructor
 */
function TShapeCalculator() {
    const service = {
        clickToCalculate,
    };

    function clickToCalculate() {
        document.getElementById('calculateBtn').addEventListener('click', startCalculation);
    }

    function startCalculation() {
        // 全部输入参数
        const inputParams = createInput();
        const calService = calculateTShapeService(inputParams);
        const allPoints = calService.calculateForceAndDelta();
        const failureMode = calService.calculateFailureMode();
        createOutput(allPoints, failureMode);
    }

    function createInput() {
        // 获取输入参数
        // const m = parseFloat(document.getElementById('m').value);
        const bf = parseFloat(document.getElementById('bf').value);
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

        // 螺栓屈服应变
        const boltQuFuEpsilon = parseFloat(document.getElementById('boltQuFuEpsilon').value);
        // 螺栓峰值应变
        const boltFengZhiEpsilon = parseFloat(document.getElementById('boltFengZhiEpsilon').value);
        // 螺栓断裂应变
        const boltDuanLieEpsilon = parseFloat(document.getElementById('boltDuanLieEpsilon').value);

        // 螺栓屈服强度
        const boltQuFuForce = parseFloat(document.getElementById('boltQuFuForce').value);
        // 螺栓峰值强度
        const boltFengZhiForce = parseFloat(document.getElementById('boltFengZhiForce').value);
        // 螺栓断裂强度
        const boltDuanLieForce = parseFloat(document.getElementById('boltDuanLieForce').value);


        // m 依赖于 翼缘厚度 和 n 值
        const m = calculateM(bf, tf, n);

        // 屈服应变
        const epsilon_y = fy / E;

        // 集中全部输入参数
        const inputParams = {
            bf, n, tf, m, lf,
            fy, E, Eh, Enk,
            boltDiameter, boltLength,
            boltHeadDiameter, washerDiameter,
            epsilon_y, epsilon_h, epsilon_m, epsilon_u,
            D_bolt, p_bolt, D_flange, p_flange,
            boltQuFuEpsilon, boltFengZhiEpsilon, boltDuanLieEpsilon,
            boltQuFuForce, boltFengZhiForce, boltDuanLieForce
        };

        // bf 翼缘总长度， 默认180mm
        // tf 翼缘厚度
        // m值 =(180-D5)/2-D3-0.8*D5
        // 依赖于 翼缘厚度 和 n 值
        function calculateM(bf, tf, n) {
            // (0.8 * tf)： 认为翼缘厚度 约等于 腹板厚度
            return (bf - tf) / 2 - n - (0.8 * tf);
        }


        return inputParams;


    }


    function calculateTShapeService(inputParams = {}) {
        const service = {
            calculateForceAndDelta, // 计算7个点的受力F及其位移Delta
            calculateFailureMode, // 计算失效模式
        };

        const {
            bf, n, tf, m, lf,
            fy, E, Eh, Enk,
            boltDiameter, boltLength,
            boltHeadDiameter, washerDiameter,
            epsilon_y, epsilon_h, epsilon_m, epsilon_u,
            D_bolt, p_bolt, D_flange, p_flange,
            boltQuFuEpsilon, boltFengZhiEpsilon, boltDuanLieEpsilon,
            boltQuFuForce, boltFengZhiForce, boltDuanLieForce
        } = inputParams;

        const  flangeService = createFlangeService();
        const  boltService = createBoltService();


        // 计算中间值
        const intermediateData = {};

        // 先算翼缘相关数据， 螺栓数据依赖翼缘的数据
        const flangeData = {
            flangeKeyCurvatures: flangeService.keyCurvatures(),
            flangeKeyMas: flangeService.calculateKeyMas(),
        }
        Object.assign(intermediateData, flangeData);

        // 螺栓相关数据， 螺栓数据依赖翼缘的数据
        Object.assign(intermediateData, {
            boltKeyPointsData: boltService.boltKeyPointsData()
        })

        const Xy = intermediateData.flangeKeyCurvatures.Xy;// 屈服曲率
        const Xh = intermediateData.flangeKeyCurvatures.Xh;// 屈服曲率
        const Xm = intermediateData.flangeKeyCurvatures.Xm;// 峰值曲率
        const Xf = intermediateData.flangeKeyCurvatures.Xf;// 断裂曲率

        const My = intermediateData.flangeKeyMas.My;// 屈服弯矩
        const Mh = intermediateData.flangeKeyMas.Mh;// 强化弯矩
        const Mm = intermediateData.flangeKeyMas.Mm;// 峰值弯矩
        const Mu = intermediateData.flangeKeyMas.Mu;// 断裂弯矩
        const Mp = intermediateData.flangeKeyMas.Mp;// 塑性弯矩


        const By = intermediateData.boltKeyPointsData.By;// 螺栓屈服荷载
        const Bu = intermediateData.boltKeyPointsData.Bu;// 螺栓峰值荷载
        const Bf = intermediateData.boltKeyPointsData.Bf;// 螺栓断裂荷载

        const Dy = intermediateData.boltKeyPointsData.Dy;// 螺栓屈服位移
        const Du = intermediateData.boltKeyPointsData.Du;// 螺栓峰值位移
        const Df = intermediateData.boltKeyPointsData.Df;// 螺栓断裂位移


        // 计算中间值相关数据 J 列（Ca），依赖螺栓相关数据
        Object.assign(intermediateData, {
            CaValues: calculateCas(), // J 列的值
        })

        // 计算中间值相关数据 H  列（X）， 曲率数据， 依赖J列数据
        Object.assign(intermediateData, {
            XValues: calculateCurvatures(), // H 列的值, 曲率列表数据
        })


        // 螺栓点的受力
        const FBy = intermediateData.CaValues.luoshanForceToDelta2.force.FBy;
        const FBu = intermediateData.CaValues.luoshanForceToDelta2.force.FBu;
        const FBf = intermediateData.CaValues.luoshanForceToDelta2.force.FBf;


        const K49 = intermediateData.CaValues.K49;
        const K50 = intermediateData.CaValues.K50;
        const J53 = intermediateData.CaValues.J53;

        // 翼缘计算service
        function createFlangeService() {
            const service = {
                keyCurvatures,
                calculateKeyMas,
                calculateTheta,
            };


// 开关函数
            function H(a, b) {
                let rs = 0;
                if ((a - b) > 0) {
                    rs = 1;
                }
                return rs;
            }

// 几个关键点的曲率
            function keyCurvatures() {

                // 屈服曲率
                const Xy = 2 * epsilon_y / tf;
                // 强化曲率
                const Xh = 2 * epsilon_h / tf;
                // 峰值曲率
                const Xm = 2 * epsilon_m / tf;
                // 断裂曲率
                const Xf = 2 * epsilon_u / tf;

                return {Xy, Xh, Xm, Xf};
            }


// 计算关键点弯矩 My,Mh, Mm, Mu,Mp
            function calculateKeyMas() {

                const {
                    bf, n, tf, m, lf,
                    fy, E, Eh, Enk,
                } = inputParams;

                const {Xy, Xh, Xm, Xf} = keyCurvatures();

                // 屈服弯矩, lf 宽度， tf厚度， fy(屈服应力410) = E(弹性模量) * epsilon_y(屈服应变)
                const My = lf * tf * tf * fy / 6;

                // 强化弯矩
                const Mh = calculateBendingMoment(Xh);

                // 峰值弯矩
                const Mm = calculateBendingMoment(Xm);

                // 断裂弯矩
                const Mu = calculateBendingMoment(Xf);

                // 塑性弯矩 (Excel D31单元格)
                const Mp = My * 3 / 2;

                return {My, Mh, Mm, Mu, Mp};
            }

// x: 点的曲率
// 返回 弯矩值
// 根据曲率计算弯矩
            function calculateBendingMoment(curvature) {

                const {
                    bf, n, tf, m, lf,
                    fy, E, Eh, Enk
                } = inputParams;

                // 屈服曲率
                // const Xy = 2 * epsilon_y / tf;
                // // 强化曲率
                // const Xh = 2 * epsilon_h / tf;
                // // 峰值曲率
                // const Xm = 2 * epsilon_m / tf;
                // // 断裂曲率
                // const Xf = 2 * epsilon_u / tf;

                const {Xy, Xh, Xm, Xf} = keyCurvatures();


                // 曲率比
                const ratioC2y = curvature / Xy;
                const ratioH2c = Xh / curvature;
                const ratioH2y = Xh / Xy;
                const ratioM2y = Xm / Xy;
                const ratioM2c = Xm / curvature;


                const item1 = ratioC2y + (0.5 * (3 - 2 * ratioC2y - Math.pow(Xy / curvature, 2)) * H(curvature, Xy));

                const item2 = 0.5 * (Eh / E) * (ratioC2y - ratioH2y) * (1 - ratioH2c) * (2 + ratioH2c) * H(curvature, Xh);

                const item3 = 0.5 * ((Eh - Enk) / E) * (ratioC2y - ratioM2y) * (1 - ratioM2c) * (2 + ratioM2c) * H(curvature, Xm);


                // 屈服弯矩, lf 宽度， tf厚度， fy(屈服应力410) = E(弹性模量) * epsilon_y(屈服应变)
                const My = lf * tf * tf * fy / 6;

                const Ma = (item1 + item2 - item3) * My;

                return Ma;


            }

            // 计算转角θ
            // Xa： 该点a的曲率X
            function calculateTheta(Xa) {
                // 计算θ的参数
                function  thetaCa() {
                    const item1 = Math.pow(Xa, 3) / (2 * Xa * Xy) -  Math.pow(Xa - Xy, 3) / (2 * Xa * Xy) * H(Xa, Xy);
                    const item2 = Eh * Math.pow(Xa - Xh, 3) / (2 * E * Xa * Xy) * H(Xa, Xh);
                    const item3 = (Eh - Enk) * Math.pow(Xa - Xm, 3) / (2 * E * Xa * Xy) * H(Xa, Xm);
                    return item1 + item2 - item3;
                }

                // ψ: J51, thy
                const thy = intermediateData.CaValues.J51;

                const Ma = calculateBendingMoment(Xa);

                const theta = m / (1 + thy) * (Xa - My / Ma * thetaCa() - 0.5 * Xy * Ma / My);

                return {theta, Ma};
            }

            return service;

        }



        // 螺栓计算service
        function createBoltService() {
            const service = {
                boltKeyPointsData,
            };

            // 计算螺栓相关参数 By, Bu, Bf
            function boltKeyPointsData() {
                const boltArea = Math.PI * boltDiameter * boltDiameter / 4;
                // 计算螺栓屈服荷载
                const By = boltArea * boltQuFuForce
                // 螺栓峰值荷载, BU(D45) = =D35^2*PI()/4*C63
                const Bu = boltArea * boltFengZhiForce;
                // 计算螺栓断裂荷载
                const Bf = boltArea * boltDuanLieForce

                // Delta y, 位移
                const Dy = boltLength * boltQuFuEpsilon; // 螺栓本身位移
                const Du = boltLength * boltFengZhiEpsilon; // 螺栓本身位移

                const boltLengthDiameterRatio = 5; // 25/5
                // 螺栓本身位移
                const Df = (boltDuanLieEpsilon * boltDiameter * boltLengthDiameterRatio)
                    - (boltFengZhiEpsilon * (boltDiameter * boltLengthDiameterRatio - boltLength));

                return  {
                    By, Bu, Bf, Dy, Du, Df
                };

            }

            return service;
        }

        // J 列，计算中间参数
        function calculateCas() {

            const D2 = m;

            // 第一步：计算基本几何参数
            const D4 = n / m;  // λ = n/m

            const D5 = tf;
            const D6 = lf;
            const D14 = E;

            const D18 = p_flange; // 翼缘率强化指数 (p)
            const D17 = D_flange;  // 翼缘率强化指数 (D)

            const D36 = 15.2; // 根据excel里面的数据
            // const D36 = 22.4; // 根据excel里面的数据
            const D38 = D_bolt || 1300000; // 螺栓率强化参数 D
            const D39 = p_bolt || 3.6; // 螺栓率强化参数 p


            // 加载速率
            const J39 = 0;
            // J40 - 需要根据Excel确定，这里假设一个值, =1+(J39/D36/D38)^(1/D39)
            const J40 = 1.0;
            const J43 = 2 * 0.1 * D2
            const J41 = J39 / 2 * D5 / J43 / D2;
            // J42 = DM (弯矩放大系数) - 需要计算
            const J42 = 1 + 2 * D18 / (1 + 2 * D18) * Math.pow((J41 / D17), 1 / D18);


            // J46 --- ξ;
            // J47 --- β;
            // J48 --- β1;
            // J49 --- β2;

            // const J46 = Mp / Mm;
            // J46: ξ = Mp/Mu，
            const J46 = intermediateData.flangeKeyMas.Mp / intermediateData.flangeKeyMas.Mu;

            // J47 = β = 2Mp/(m*Bu)
            const J47 = 2 * intermediateData.flangeKeyMas.Mp / (m * intermediateData.boltKeyPointsData.Bu);

            // J48, J49 - 失效模式阈值
            // J48 --- β1; J49 --- β2;
            const J48 = 2 * D4 * J46 / ((1 + 2 * D4) * J42);;
            const J49 = calculateJ49Threshold(J40, J46, J42, D4);

            // 第五步：计算失效模式J50
            const J50 = calculateJ50(J47, J48, J49, J40);

            //  计算失效模式过渡系数 ψ
            const J51 = calculateJ51(J47, D4, J46, J42, J40);

            const J52 = calculateJ52(D4, J51);


            // 第八步：最终计算K49和K50
            const K49 = calculateK49(J47, J48, J49);
            // const K50 = calculateK50(J50, J47, J48, J49, N49, m, tf, D42);
            const K50 = calculateK50(J47, J40);

            // 螺栓刚度：(螺栓材料的强度)
            const J53 = calculateJ53(D14, D6, D5, D2);

            // 翼缘速率 mm/s
            const J54 = calculateJ54(J39, J47, J49, K49, D4);

            //
            const J56 = calculateJ56(J39, J47, J49, K49, D4);

            // 翼缘边缘应变率
            const J57 = calculateJ57(J54, D5, J43, D2);

            // 率强化修正系数
            const J58 = calculateJ58(D18, J57, D17);

            // J59
            const J59 = calculateJ59(J56, D36, D38, D39);

            // 螺栓 F - Δ2 （螺栓受力位移）, 此数据用作插值参考
            const luoshanForceToDelta2 = calculateBoltForceAndMove(J52, J59);


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

            // J52计算函数
            function calculateJ52(D4, J51) {
                // Excel公式: =2*D4*(1+J51)/(D4*(1+J51)+J51)

                const numerator = 2 * D4 * (1 + J51);
                const denominator = D4 * (1 + J51) + J51;

                if (denominator === 0) {
                    return 2; // 默认值
                }

                return numerator / denominator;
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

            // J56计算函数
            // J56=J39*IF(J47>J49,1-K49/(1+D4),D4*K49/(1+D4))
            function calculateJ56(J39, J47, J49, K49, D4) {
                // 参数验证
                if (typeof J39 !== 'number' || typeof J47 !== 'number' ||
                    typeof J49 !== 'number' || typeof K49 !== 'number' ||
                    typeof D4 !== 'number') {
                    // console.error("J56计算错误: 所有参数必须为数字");
                    return 0;
                }

                // 避免除零错误
                if (1 + D4 === 0) {
                    // console.warn("J56计算警告: 1+D4为0，返回0");
                    return 0;
                }

                let conditionTerm;

                if (J47 > J49) {
                    // 情况1: J47 > J49
                    conditionTerm = 1 - K49 / (1 + D4);
                    // console.log(`J56计算: J47(${J47}) > J49(${J49}), 使用公式: 1 - ${K49}/(1+${D4}) = ${conditionTerm}`);
                } else {
                    // 情况2: J47 <= J49
                    conditionTerm = (D4 * K49) / (1 + D4);
                    // console.log(`J56计算: J47(${J47}) <= J49(${J49}), 使用公式: ${D4}*${K49}/(1+${D4}) = ${conditionTerm}`);
                }

                const result = J39 * conditionTerm;
                // console.log(`J56最终结果: ${J39} * ${conditionTerm} = ${result}`);

                return result;
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
                if (base < 0 && 1 / D18 % 2 === 0) {
                    return 1; // 如果底数为负且指数分母为偶数，返回基准值
                }

                // 计算幂
                const powerTerm = Math.pow(Math.abs(base), 1 / D18);

                // 如果底数为负且指数分母为奇数，处理符号
                let signedPowerTerm = base < 0 ? -powerTerm : powerTerm;

                // 最终计算
                return 1 + coefficient * signedPowerTerm;
            }

            // J59计算函数
            function calculateJ59(J56, D36, D38, D39) {
                // Excel公式: =1+(J56/D36/D38)^(1/D39)

                // 避免除零错误
                if (D36 === 0 || D38 === 0 || D39 === 0) {
                    return 1;
                }

                const normalizedRate = J56 / D36 / D38;

                if (normalizedRate < 0) {
                    return 1;
                }

                const rateEffect = Math.pow(normalizedRate, 1 / D39);
                return 1 + rateEffect;
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


            // 计算螺栓载荷及其位移数据
            function calculateBoltForceAndMove(J52, J59) {

                // 从螺栓屈服荷载算出力F
                const FBy = By * J59 * J52;

                console.log(`
                    By * J59 * J52: ${By} * ${J59} * ${J59} 
                `);

                // 从螺栓峰值荷载算出力F
                const FBu = Bu * J59 * J52;

                // 从螺栓断裂荷载算出力F
                const FBf = Bf * J59 * J52;


                return {
                    force: {FBy, FBu, FBf},
                    move: {Dy, Du, Df},
                };
            }


            return {
                J39, J40,  J41, J42, J43, J46, J47, J48, J49, J50,  J51 , J52, J53,  J54 , J56, J57,  J58 , J59,
                K49, K50,
                luoshanForceToDelta2,
            };
        }


        // H 列， 计算曲率列的值及其插值
        function calculateCurvatures () {

            // 在两个元素之间均等插入 n 个值
            // 在arr这个数组的 start和end之间等距插入n个值
            function insertBetween(arr, start, end, n) {
                if (!arr.includes(start) || !arr.includes(end)) {
                    throw new Error('start 或 end 不在数组中');
                }

                const diff = end - start;
                const step = diff / (n + 1);

                const values = [];
                for (let i = 1; i <= n; i++) {
                    values.push(start + step * i);
                }

                const result = [...arr];
                const endIndex = result.indexOf(end);
                result.splice(endIndex, 0, ...values);

                return result;
            }

            // 曲率数组列表，从Xy屈服曲率开始
            let Xs = [ Xy, Xh,  Xm, Xf, ];
            // h 和 m之间插入4个值
            Xs = insertBetween(Xs, Xh, Xm,4);
            // m 和 f之间插入19个值
            Xs = insertBetween(Xs, Xm, Xf,19);

            // K 列， 转角值列表
            const thetas = [];

            // L 列， 受力插值列表
            const predictForces = [];
            // ψ: J51, thy
            const thy = intermediateData.CaValues.J51;
            Xs.map(X => {
                const {theta, Ma} = flangeService.calculateTheta(X);
                thetas.push(theta);
                const force = 2 * Ma * (1 + thy) / ( m * Math.cos(theta));
                predictForces.push(force);
            });

            return {Xs, thetas, predictForces};
        }


        // 计算 T - 实现TREND函数逻辑
        // R(F): 为受力（R列的值）， 参考插值表的值，用trend（excel函数）插值算法根据Ma（弯矩的值）得到 T（Δ2）的值
        // Δ2， T 列， 螺栓产生的位移
        function flangeDelta2(Rf) {
            // const luoshanForceToDelta2 = intermediateData.CaValues.luoshanForceToDelta2;
            // 前面的 + 号把 toFixed 后的字符串转换回来成数字
            // const FBy = +FBy.toFixed(0);
            // const FBu = +FBu.toFixed(0);
            // const Dy = +Dy.toFixed(2);
            // const Du = +Du.toFixed(2);

            const xs = [0, FBy, FBu];
            const ys = [0, Dy, Du];

            // console.log(xs, ys)

            function chazhi(x, x1, y1, x2, y2) {
                let y = (x - x1) * (y2 - y1) / (x2 - x1) + y1
                return y;
            }

            let result = 0;
            if (Rf <= xs[1]) {
                result = chazhi(Rf, xs[0], ys[0], xs[1], ys[1]);
            } else {
                result = chazhi(Rf, xs[1], ys[1], xs[2], ys[2]);
            }
            return result;
        }

        // Δ1， S 列， 翼缘产生的位移
        // 依赖 delta2， Δ2 （T列）
        // R(f): 为受力（R列的值）
        function flangeDelta1(Rf, theta, delta2) {


            // λ: = n/m
            const D4 = n / m;

            // S5=$K$49*SIN(K5)*$D$2+($K$50-$K$49)*T5/$D$4/2+R5/$J$53
            // const S = K49 * Math.sin(theta) * D2 + (K50 - K49) * theta / D4 / 2 + R / J53;
            const S = K49 * Math.sin(theta) * m + (K50 - K49) * delta2 / D4 / 2 + Rf / J53;
            return S;
        }

        // 螺栓产生的Δ1位移（翼缘）
        // R(f): 为受力（R列的值）
        // 依赖 delta2， Δ2 （T列）
        function boltDelta1 (Rf, delta2) {

            // Xs: 曲率, thetas： 转角, predictForces： 受力
            const {Xs, thetas, predictForces} = calculateCurvatures();
            // console.log({Xs, thetas, predictForces});

            // 通过插值预测出来的 转角 θ (theta)
            // 从力到转角
            const boltPredictTheta = excelTrend(thetas, predictForces, Rf);

            // λ: = n/m
            const D4 = n / m;
            // 螺栓产生的Δ1位移, S列
            const boltDelta1 = K49 * Math.sin(boltPredictTheta) * m + (K50 - K49) * delta2 / D4 / 2 + Rf / J53;

            return boltDelta1;

        }


        // 计算螺栓断裂位移
        // S9=S8+(R9-R8)/$J$53
        // Su（S8）:  螺栓峰值(Sf 断裂的前一个状态)时引起的翼缘位移 S
        // Ff(R9):  断裂状态的受力F
        // Fu(R8):  螺栓峰值(Sf 断裂的前一个状态)时 的受力
        // T，螺栓本身的位移，T列的值
        function boltBfDelta1_bak(Su, Ff, Fu, T) {
            const J53 = intermediateData.CaValues.J53;
            const S = Su + (Ff - Fu) / J53;
            const U = 2 * S + T;
            return U;
        }
        function boltBfForceMove(BuDelta1, FBf, FBu, BfDelta2) {
            const J53 = intermediateData.CaValues.J53;
            const d2 = BfDelta2;
            const force = FBf;
            const d1 = BuDelta1 + (FBf - FBu) / J53;
            const move = 2 * d1 + d2;

            return {force, move, d1, d2};
        }

        // 计算7个点的受力F及其位移Delta
        function calculateForceAndDelta() {

            // 计算 R (翼缘  F)
            function flangeForceAndMove(moment,  theta) {
                // ψ: J51, thy
                const thy = intermediateData.CaValues.J51;
                const force = 2 * moment * (1 + thy) / ( m * Math.cos(theta));

                // 螺栓位移 delta2 (T)
                const d2 = flangeDelta2(force);

                // 翼缘位移 delta1 (S)
                const d1 = flangeDelta1(force, theta, d2);

                const move = 2 * d1 + d2;

                return {force, move};
            }

            function boltForceAndMove(boltPoint) {
                // 螺栓位移 delta2 (T)
                const d2 = boltPoint.delta2;
                const force = boltPoint.force;
                const d1 = boltDelta1 (force, d2);
                const move = 2 * d1 + d2;

                return {force, move, d1, d2};
            }


            // 翼缘点的处理
            // c: 曲率, Ma: 弯矩
            let flangePoints = [
                { x: 0, y: 0, name: "翼缘屈服点", id: 'May', c: Xy, Ma: My },
                { x: 0, y: 0, name: "翼缘强化点", id: 'Mah', c: Xh, Ma: Mh },
                { x: 0, y: 0, name: "翼缘峰值点", id: 'Mam', c: Xm, Ma: Mm },
                { x: 0, y: 0, name: "翼缘断裂点", id: 'Mau', c: Xf, Ma: Mu },
            ];
            flangePoints = flangePoints.map(i => {
                const {theta, Ma} = flangeService.calculateTheta(i.c);
                const forceMove = flangeForceAndMove(i.Ma, theta);
                i.x = forceMove.move;
                i.y = forceMove.force;
                return i;
            });


            // FB: 载荷
            let boltPoints = [
                { x: 0, y: 0, name: "螺栓屈服点", id: 'By', force: FBy, delta2: Dy },
                { x: 0, y: 0, name: "螺栓峰值点", id: 'Bu', force: FBu, delta2: Du },
                { x: 0, y: 0, name: "螺栓断裂点", id: 'Bf', force: FBf, delta2: Df },
            ];

            const assignBoltValues = (point, forceMove) => {
                Object.assign(point, {
                    x: forceMove.move,
                    y: forceMove.force,
                    d1: forceMove.d1,
                    d2: forceMove.d2,
                });
            };
            // 计算三个点的 force和位移
            assignBoltValues(boltPoints[0], boltForceAndMove(boltPoints[0])); // Bu
            assignBoltValues(boltPoints[1], boltForceAndMove(boltPoints[1])); // By
            assignBoltValues(boltPoints[2], boltBfForceMove(boltPoints[1].d1, FBf, FBu, Df)); // Bf


            console.log(
            'boltPoints: ' ,boltPoints,
            'flangePoints： ', flangePoints,
            );


            const the7points = sortThe7Points(boltPoints, flangePoints);

            return the7points;

        }


        /**
         * @desc  Excel TREND函数的JavaScript实现
         *
         * knownYs: Y值数组 [y1, y2, ...]
         * knownXs: X值数组 [x1, x2, ...]
         * newX: 要插值的X值
         * 返回: 对应的Y值
         */
        function excelTrend(knownYs, knownXs, newX) {


            // 参数验证
            if (!Array.isArray(knownYs) || !Array.isArray(knownXs)) {
                throw new Error("knownYs和knownXs必须是数组");
            }

            if (knownYs.length !== knownXs.length) {
                throw new Error("knownYs和knownXs数组长度必须相同");
            }

            if (knownYs.length < 2) {
                throw new Error("至少需要2个数据点进行插值");
            }

            // 如果newX正好等于某个已知点，直接返回对应的Y值
            for (let i = 0; i < knownXs.length; i++) {
                if (knownXs[i] === newX) {
                    return knownYs[i];
                }
            }

            // 找到newX所在的区间
            let lowerIndex = -1;
            let upperIndex = -1;

            // 对已知点按X值排序（保持X-Y对应关系）
            const points = knownXs.map((x, index) => ({x, y: knownYs[index]}));

            // 不排序了，就按参考的数据
            points.sort((a, b) => a.x - b.x);

            const sortedXs = points.map(p => p.x);
            const sortedYs = points.map(p => p.y);

            // 处理边界情况
            if (newX <= sortedXs[0]) {
                // newX小于等于最小值，使用第一个区间
                lowerIndex = 0;
                upperIndex = 1;
            } else if (newX >= sortedXs[sortedXs.length - 1]) {
                // newX大于等于最大值，使用最后一个区间
                lowerIndex = sortedXs.length - 2;
                upperIndex = sortedXs.length - 1;
            } else {
                // 在中间找到newX所在的区间
                for (let i = 0; i < sortedXs.length - 1; i++) {
                    if (newX >= sortedXs[i] && newX <= sortedXs[i + 1]) {
                        lowerIndex = i;
                        upperIndex = i + 1;
                        break;
                    }
                }
            }

            // 如果没找到合适的区间（理论上不会发生）
            if (lowerIndex === -1 || upperIndex === -1) {
                throw new Error("无法找到合适的插值区间");
            }

            // 线性插值计算
            const x1 = sortedXs[lowerIndex];
            const y1 = sortedYs[lowerIndex];
            const x2 = sortedXs[upperIndex];
            const y2 = sortedYs[upperIndex];


            // 避免除零
            if (x1 === x2) {
                return (y1 + y2) / 2; // 如果X值相同，返回Y值的平均值
            }

            // 线性插值公式: y = y1 + (newX - x1) * (y2 - y1) / (x2 - x1)
            const result = y1 + (newX - x1) * (y2 - y1) / (x2 - x1);

            // if([37779, 78255].includes(parseInt(newX))) {
            //     // console.log("37779 ==== ",, x2, y2, newX)
            //     console.log(`
            //      插值：
            //         newX: ${newX}
            //         (x1, y1): ${x1}, ${y1},
            //         (x2, y2): ${x2}, ${y2},
            //
            //         result: ${result}
            //     `);
            // }

            return result;
        }



        /**
         *
         *  对7个点[除了螺栓的最后一个点(Δ，Fbf)]的纵坐标从小到大排序，
         *  如果从小到大，若该点纵坐标达到Fu，后续数据就不再取了；
         *  若该点数据达到Fbu，则下一点为Fbf，若该点为Fbf，后续数据就不再取了
         *
         *  Fu(翼缘弯矩断裂) 理论上应该少于Fbu(螺栓峰值载荷)，按从小到大的排序，到了Fbu点之后的下一个肯定是Fbf点（螺栓断裂点载荷）, 都断裂了，展示该点没意义
         *
         *  ... Fu ... Fbu Fbf ...  (Fu < Fbu)
         *  ... Fbu Fbf ... Fu ...  (Fu > Fbu)
         *
         * 这里边是7个点，但是最后排列的时候有一定的规则[捂脸]两种情况：
         * 1、从小到大排列（不包括Bf），超过Mu对应的荷载的点，都不再参与绘制曲线；
         * 2、从小到大排列（不包括Bf），超过Bu的点都不再参与绘制曲线，Bu后一个点必须是Bf
         *
         *  @var boltPoints , 螺栓计算出来的3个关键点，按顺序是 （δy，By）、（δu，Bu）、（δf，Bf）
         *  @var flangePoints , 翼缘计算出来的4个关键点，按顺序是 （θy，May）、（θh，Mah）、（θm，Mam）、（θu，Mau）
         *
         *
         *  { x: delta1, y: F1, name: "螺栓屈服点", id: 'By' },
         *  { x: delta2, y: F2, name: "螺栓峰值点", id: 'Bu'  },
         *  { x: delta3, y: F3, name: "螺栓断裂点", id: 'Bf'  }
         *
         *  { x: thePoint.U, y: thePoint.R  / 1000, name: "翼缘屈服点", id: 'May' }
         *  { x: thePoint.U, y: thePoint.R  / 1000, name: "翼缘强化点", id: 'Mah' }
         *  { x: thePoint.U, y: thePoint.R  / 1000, name: "翼缘峰值点", id: 'Mam' }
         *  { x: thePoint.U, y: thePoint.R  / 1000, name: "翼缘断裂点", id: 'Mau' }
         *
         *
         */
        function sortThe7Points(boltPoints, flangePoints) {


            const point00 = {x: 0, y: 0, name: '原点', id: 'origin'};

            // 合并所有点并按纵坐标排序, 去掉Bf点
            const allPoints = [point00, ...flangePoints, ...boltPoints].sort((a, b) => a.y - b.y).filter(i => i.id !== 'Bf');

            // Fu(翼缘弯矩断裂点)
            const pointFu = allPoints.find(i => i.id === 'Mau') || {};

            // Fbu(螺栓峰值载荷点,其下一点必定是螺栓断裂点Bf)
            const pointFbu = allPoints.find(i => i.id === 'Bu') || {};

            // Fbf(螺栓断裂点Bf)
            const pointFbf = boltPoints.find(i => i.id === 'Bf') || {};


            let drawPoints = [];
            // Fu < Fbu 的情况，点只截取到Fu为止
            if (pointFu.y <= pointFbu.y) {

                // 截取到第一个弯矩断裂点为止
                const pointFuIndex = allPoints.findIndex(i => i.id === 'Mau');
                drawPoints = allPoints.slice(0, pointFuIndex + 1);

            } else {
                // Fu > Fbu 的情况，点只截取到Fbu为止

                // 截取到第一个螺栓峰值点为止（其下一个点就是断裂点）
                const pointFbuIndex = allPoints.findIndex(i => i.id === 'Bu');
                drawPoints = allPoints.slice(0, pointFbuIndex + 1);

                //因为最后的一个点是Bu, 顺带把Bf点也带上
                drawPoints.push(pointFbf);
            }

            return drawPoints;

        }


        // 计算失效模式
        function calculateFailureMode() {
            const failureMode = intermediateData.CaValues.J50 || 'FM1-FF';
            const J47 = intermediateData.CaValues.J47.toFixed(3);
            const J51 = intermediateData.CaValues.J51.toFixed(3);
            const params = `<span class="f-model-b">β = ${J47}</span>,  <span class="f-model-b">ψ = ${J51}</span> `;
            return `<span class="f-model">${failureMode}</span>, ${params}`;
        }

        return service;

    }

    function createOutput(allPoints, failureMode) {

        allPoints = useKnUnit();

        // 页面加载时初始化图表
        window.addEventListener('DOMContentLoaded', function () {
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

        // 显示结果
        displayResults(allPoints, failureMode);

        // 绘制图表
        // drawChart(allPoints, failureMode);
        drawChart(allPoints, failureMode);

        // 显示计算结果
        function displayResults(points, failureMode) {
            // 显示失效模式
            document.getElementById('failureModeResult').innerHTML = `<strong>${failureMode}</strong>`;

            // 显示关键数据点
            let pointsHTML = '';
            points.forEach((point, index) => {
                // 原点虽然画出来，但是不标记，下边也不展示关键信息
                if (index > 0) {
                    pointsHTML += `
                <div class="result-item">
                     <div class="r-i1"><strong>${index}, ${point.name}  ${point.id} </strong>:</div>
                     <div class="r-i1">荷载F = ${point.originalY.toFixed(0)} N</div>
                     <div class="r-i1">变形Δ = ${point.x.toFixed(2)} mm</div>
                </div>
            `;
                }
            });
            document.getElementById('pointsResult').innerHTML = pointsHTML;
        }

        // 用千牛做单位
        // 返回全部点
        function useKnUnit() {
            return allPoints.map(i => {
                i.originalY = i.y;
                i.y = i.y / 1000;
                return i;
            });
        }


        function drawChart(points) {
            const chartDom = document.getElementById('chart');
            const myChart = echarts.init(chartDom);

            const option = {
                title: {
                    text: 'T形件全过程曲线',
                    left: 'center'
                },
                tooltip: {
                    trigger: 'item',
                    formatter: function (params) {
                        if (params.dataIndex === 0) {
                            return ''; // 原点不展示tooltip
                        }
                        const point = points[params.dataIndex];
                        return `
                    ${point.name} ${point.id}<br/>
                    荷载F: ${point.y.toFixed(2)} kN <br> 
                    变形Δ: ${point.x.toFixed(4)} mm
                `;
                    }
                },
                xAxis: {
                    type: 'value',
                    name: '变形 Δ (mm)'
                },
                yAxis: {
                    type: 'value',
                    name: '荷载 F (kN)'
                },
                series: [{
                    type: 'line',
                    data: points.map(point => [point.x, point.y]),
                    // showSymbol: 'circle',
                    // showSymbol: false,
                    // symbolSize: 8,
                    lineStyle: {
                        color: '#3498db',
                        width: 2
                    },

                    markPoint: {

                        data: points.map((point, index) => {
                            let mPoint = {
                                name: `${index + 1}`,
                                //coord: [point.x , point.y],
                                coord: [point.x, point.y + 1.5],
                                symbol: 'circle',
                                symbolSize: 5,
                                itemStyle: {
                                    color: '#fff',
                                    borderColor: '#e74c3c',
                                    borderWidth: 0
                                },
                                label: {
                                    show: true,
                                    // formatter: '{b}',
                                    // 原点虽然画出来，但是不标记，下边也不展示关键信息
                                    formatter: function (params) {
                                        // console.log('formatter params:', params);

                                        const pointIndex = params.dataIndex;
                                        // const point = filteredPoints[pointIndex];

                                        // 1. 跳过特定点的数字标记（如果需要）
                                        if (pointIndex === 0) { // 原点
                                            return ''; // 返回空字符串则不显示数字
                                        } else {
                                            return pointIndex;
                                        }

                                    },
                                    position: 'inside',
                                    color: '#000',
                                    fontWeight: 'normal'
                                }
                            }
                            if (index === 0) {
                                // 原点不展示
                                Object.assign(mPoint, {
                                    symbol: '',
                                    symbolSize: 0,
                                    label: {
                                        show: false
                                    }
                                });
                            }

                            return mPoint;
                        })


                    }

                }],
                grid: {
                    left: '10%',
                    right: '5%',
                    bottom: '15%',
                    top: '15%'
                }
            };

            myChart.setOption(option);

            window.addEventListener('resize', function () {
                myChart.resize();
            });
        }
    }

    return service;
}

const tCalculator = TShapeCalculator();
tCalculator.clickToCalculate();