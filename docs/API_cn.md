# coligunCalc C++ API 参考

> English version: [API.md](API.md)

---

## 快速开始

### 引入与链接

整个库位于 `namespace coilgun` 下。最简单的引入方式：

```cpp
#include <coilgun/coilgun.hpp>   // 全部物理层 + 组件层 + 仿真层

// 或者按需引入：
#include <coilgun/physics/constants.hpp>
#include <coilgun/components/driving_coil.hpp>
#include <coilgun/simulation/single_stage_sim.hpp>
// ...
```

链接静态库：

```cmake
target_link_libraries(your_target PRIVATE coilgun)
```

临时脚本也可以直接编译链接：

```sh
g++ -std=c++17 -Iinclude your_file.cpp build/src/libcoilgun.a -o your_binary
```

### 命名空间

| 命名空间 | 内容 |
|-----------|----------|
| `coilgun::physics` | 物理常量、椭圆积分、Struve 函数、数值求积、自感/互感、LRU 缓存、查表 |
| `coilgun::components` | DrivingCoil 和 Armature 类 |
| `coilgun::simulation` | 仿真引擎：时间步进器、激励模型、终止策略、触发配置、SimState/MultiStageState、SingleStageSim、MultiStageSim |
| `coilgun::physics::detail` | 内部实现细节（查表数据）——不应直接依赖 |

### 最小示例

计算铜质驱动线圈的自感和直流电阻：

```cpp
#include <coilgun/coilgun.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;

    // 150匝铜线圈：内径10mm、外径30mm、长度50mm、线径1mm²、填充率70%
    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);

    std::cout << "自感：  " << coil.self_inductance() << " H\n";
    std::cout << "直流电阻：" << coil.resistance() << " Ω\n";
    std::cout << "匝密度：  " << coil.turns_density() << " 匝/m²\n";
    return 0;
}
```

计算两个同轴电流丝环之间的互感：

```cpp
#include <coilgun/physics/mutual_inductance.hpp>
#include <iostream>

int main() {
    double M  = coilgun::physics::mutual_inductance_filament(0.02, 0.025, 0.01);
    double dM = coilgun::physics::mutual_inductance_gradient_filament(0.02, 0.025, 0.01);

    std::cout << "M  = " << M  << " H\n";
    std::cout << "dM = " << dM << " H/m\n";
    return 0;
}
```

---

## 架构

```
include/coilgun/
├── core/types.hpp              — 公共类型、前向声明
├── physics/
│   ├── constants.hpp           — μ₀、材料属性（铜/铝）、cp(T)、rho(T)
│   ├── elliptic.hpp            — 第一/二类完全椭圆积分 K(k)、E(k)
│   ├── struve.hpp              — Struve 函数 H₀(x)、H₁(x)
│   ├── quadrature.hpp          — Gauss-Legendre 与 Gauss-Laguerre 节点
│   ├── cache.hpp               — LRU 缓存模板（默认 4096 条）
│   ├── lookup_tables.hpp       — T(q,p) 电感形状因子查表
│   ├── lookup_table_data.hpp   — 预计算 T(q,p) 稠密表（约 290 万条）
│   ├── self_inductance.hpp     — 自感（精确计算 + 快速查表）
│   └── mutual_inductance.hpp   — 丝级和线圈级互感
├── components/
│   ├── driving_coil.hpp        — 驱动线圈类
│   └── armature.hpp            — 电枢类（m×n 电流丝离散化）
├── simulation/
│   ├── excitation.hpp          — Excitation、CapacitorExcitation、CrowbarExcitation、WaveformExcitation
│   ├── time_stepper.hpp        — EulerStepper、RK4Stepper
│   ├── sim_result.hpp          — SimStep、SimResult、SimSummary
│   ├── termination.hpp         — TerminationPolicy
│   ├── single_stage_sim.hpp    — SimState、SingleStageSim<StepperPolicy>
│   ├── trigger_config.hpp      — TriggerMode、TriggerConfig
│   ├── multi_stage_result.hpp  — StepSnapshot、MultiStageStep、PerStageSummary、MultiStageSummary、MultiStageResult
│   └── multi_stage_sim.hpp     — OptimizationLevel、MultiStageState、MultiStageSim<StepperPolicy>
└── coilgun.hpp                 — 便利总头文件
```

---

## 物理常量

```cpp
#include <coilgun/physics/constants.hpp>

namespace coilgun::physics {
    constexpr double MU0 = 4.0 * M_PI * 1e-7;    // 真空磁导率 (H/m)
    constexpr double T_REFERENCE = 293.0;          // 室温 (K)

    struct MaterialProperties {
        double resistivity_ref;    // 参考温度电阻率 (Ω·m)
        double temp_coefficient;   // 温度系数 β (K⁻¹)
        double density;            // 质量密度 (kg/m³)
    };

    extern const MaterialProperties COPPER;        // 1.75e-8 Ω·m, 0.0041 K⁻¹, 8960 kg/m³
    extern const MaterialProperties ALUMINUM;      // 2.82e-8 Ω·m, 0.0042 K⁻¹, 2700 kg/m³

    enum class ArmatureMaterial { Aluminum, Copper };

    double material_cp(ArmatureMaterial m, double T);    // 比热容调度
    double material_beta(ArmatureMaterial m);             // 电阻率温度系数调度

    double specific_heat_capacity_copper(double T);       // J/(kg·K), 式 6.3
    double specific_heat_capacity_aluminum(double T);     // J/(kg·K), 式 6.2
    double resistivity_copper(double T);                  // Ω·m, 式 6.9
    double resistivity_aluminum(double T);                // Ω·m, 式 6.9
}
```

**用法**：

```cpp
double rho_cu_at_373K = coilgun::physics::resistivity_copper(373.0);
// ≈ 2.28e-8 Ω·m（铜在 100°C 时的电阻率）

// 将材料属性传给组件：
DrivingCoil coil(..., COPPER.resistivity_ref, ...);
```

`material_cp` 根据枚举值分派 `specific_heat_capacity_copper` 或 `specific_heat_capacity_aluminum`。`material_beta` 分派温度系数。

---

## 椭圆积分

```cpp
#include <coilgun/physics/elliptic.hpp>

double coilgun::physics::elliptic_k(double m);    // 第一类完全椭圆积分 K(m)
double coilgun::physics::elliptic_e(double m);    // 第二类完全椭圆积分 E(m)
double coilgun::physics::elliptic_modulus(
    double radius_a, double radius_b,
    double separation);                           // 两同轴环的椭圆模量 k
```

封装 Boost.Math 的 `ellint_1` / `ellint_2`。使用**参数约定**（m = k²），而非模量约定。

`elliptic_modulus` 计算两个同轴圆环在给定间距下的几何模量——所有丝级互感计算的基础（NumericalModel.md 式 4.7）：

```
k = sqrt( 4·a·b / ((a + b)² + h²) )
```

其中 `a`、`b` 为两个环的半径，`h` 为轴向间距。

**注意**：内部对 `k` 做了截断，避免 `k → 1` 时的奇点（重合环），因此 `elliptic_modulus(a, a, 0)` 会返回略小于 1 的有限值。

**快速验证**（已知值）：

```cpp
elliptic_k(0.0);   // π/2 ≈ 1.570796
elliptic_e(0.0);   // π/2 ≈ 1.570796
elliptic_k(0.5);   // ≈ 1.854075
elliptic_e(0.5);   // ≈ 1.350644
```

---

## Struve 函数

```cpp
#include <coilgun/physics/struve.hpp>

double coilgun::physics::struve_h0(double x);     // Struve H₀(x)
double coilgun::physics::struve_h1(double x);     // Struve H₁(x)
```

**奇偶性**：H₀ 为奇函数（H₀(-x) = -H₀(x)），H₁ 为偶函数（H₁(-x) = H₁(x)）。

SciPy 三区间策略，按 |x| 分段：
- **|x| < 8**：幂级数——快速，机器精度级准确
- **|x| ≥ 20**：渐进展开，通过修正 Bessel K₀/K₁ 积分 + Gauss-Laguerre 求积
- **8 ≤ |x| < 20**：两种方法同时计算，结果交叉验证至 1e-6 相对容差

已与 `scipy.special.struve` 验证。由 `self_inductance_exact` 内部调用；通常不需要直接使用。

---

## 数值求积

```cpp
#include <coilgun/physics/quadrature.hpp>

struct coilgun::physics::QuadratureNodes {
    std::vector<double> nodes;    // 求积节点位置
    std::vector<double> weights;  // 求积权重
};

QuadratureNodes gauss_legendre(int n);    // n 点 Gauss-Legendre 规则，区间 [-1, 1]
QuadratureNodes gauss_laguerre(int n);    // n 点 Gauss-Laguerre 规则，区间 [0, +∞)
```

内部采用 Golub-Welsch 算法（Jacobi 矩阵特征值）。被 `self_inductance_exact`、`mutual_inductance_coil` 和 Struve 渐进区间调用。普通用户通常不需要直接使用。

---

## 自感

```cpp
#include <coilgun/physics/self_inductance.hpp>

// 主要 API — 根据几何参数自动查表，越界时回落精确计算

double coilgun::physics::self_inductance(
    double inner_radius, double outer_radius, double length,
    double turns_density,
    bool force_exact = false);    // (H) — 默认查表，true 则跳过查表强制精确

// 当需要显式获得参考级精度时使用

double coilgun::physics::self_inductance_exact(
    double inner_radius, double outer_radius, double length,
    double turns_density);                        // (H) — 始终使用 Bessel/Struve 核 + GL16

// 形状因子（更低层，无量纲）

double coilgun::physics::inductance_shape_factor(double q, double p);          // 查表 + 双线性插值
double coilgun::physics::inductance_shape_factor_reference(double q, double p); // 精确 Bessel/Struve（慢）
```

空心圆柱线圈的自感为：

```
L = 2π · μ₀ · nc² · ri⁵ · T(q, p)
```

其中 `nc = N / ((ro - ri) × l)` 为匝密度，`ri` 为内半径，`T(q, p)` 为无量纲形状因子，`q = 长度 / 内半径`，`p = 外半径 / 内半径`。

**`self_inductance` — 推荐默认使用。**

当 `force_exact` 为 false（默认）且 `(q, p)` 落在查表范围内时（`q ∈ [0.05, 4]`，`p ∈ [1.05, 4]`），使用快速 T(q,p) 查表（~μs，T 值 ~5e-7 相对误差）。超出范围则自动回落至精确 Bessel/Struve 积分（~ms）。当 `force_exact` 为 true 时，跳过查表无条件使用参考级积分。

**`self_inductance_exact` — 参考级精度保证。**

始终使用完整的 Bessel/Struve 核 + 复合 GL16 积分。等价于调用 `self_inductance` 时设置 `force_exact = true`。

| 函数 | 行为 | 速度 | 适用场景 |
|----------|-----------|-------|----------|
| `self_inductance` | 范围内查表，范围外精确 | ~μs 或 ~ms | **默认** — 仿真、组件 |
| `self_inductance_exact` | 始终精确 | ~ms | 验证、基准数据 |

`DrivingCoil` 和 `Armature` 类内部均使用 `self_inductance`。

**示例**：

```cpp
// q = 0.05/0.01 = 5 超出查表范围，因此回落到精确计算：
double L1 = self_inductance(0.01, 0.03, 0.05, 1e5);

// q = 0.0002/0.01 = 0.02 < 0.05 — 触发精确回落：
double L2 = self_inductance(0.01, 0.03, 0.0002, 1e5);

// 显式请求参考级：
double L3 = self_inductance_exact(0.01, 0.03, 0.05, 1e5);

// 通过 bool 标志强制精确（与 self_inductance_exact 等价）：
double L4 = self_inductance(0.01, 0.03, 0.05, 1e5, true);
```

---

## 互感

```cpp
#include <coilgun/physics/mutual_inductance.hpp>
```

### 丝级——解析表达式

两个同轴圆环。这是最内层核函数；调用结果有 LRU 缓存（4096 条）。

```cpp
double mutual_inductance_filament(
    double radius_a, double radius_b, double separation);     // (H), 式 4.8

double mutual_inductance_gradient_filament(
    double radius_a, double radius_b, double separation);     // (H/m), 式 4.15
```

`separation` 为两个环平面之间的轴向距离。梯度 `dM/dz` 关于间距是**奇函数**：间距符号翻转，梯度符号随之翻转。

**示例**：两个半径 20mm 和 25mm 的环，间距 10mm：

```cpp
double M  = mutual_inductance_filament(0.02, 0.025, 0.01);   // nH 量级
double dM = mutual_inductance_gradient_filament(0.02, 0.025, 0.01);
```

### 线圈级——4D Gauss-Legendre 求积

完整有限截面的厚线圈。每维使用 n 点 Gauss-Legendre（共 n⁴ 次核函数计算），每个内部丝级调用均经过 LRU 缓存。

```cpp
double mutual_inductance_coil(
    double r_inner_a, double r_outer_a, double length_a, int turns_a,
    double r_inner_b, double r_outer_b, double length_b, int turns_b,
    double separation,
    int n_nodes = 9);                                       // (H), 式 4.13

double mutual_inductance_gradient_coil(
    double r_inner_a, double r_outer_a, double length_a, int turns_a,
    double r_inner_b, double r_outer_b, double length_b, int turns_b,
    double separation,
    int n_nodes = 9);                                       // (H/m), 式 4.16
```

此处的 `separation` 为线圈中心到中心的轴向距离（与丝级的面到面距离不同）。

**示例**——驱动线圈与模拟为厚线圈的电枢段之间的互感：

```cpp
double M_coil = mutual_inductance_coil(
    0.01, 0.03, 0.05, 150,    // 线圈 A: ri, re, len, 匝数
    0.005, 0.025, 0.04, 1,    // 线圈 B: ri, re, len, 匝数（电枢段）
    0.06);                     // 中心到中心间距

double dM_coil = mutual_inductance_gradient_coil(
    0.01, 0.03, 0.05, 150,
    0.005, 0.025, 0.04, 1,
    0.06);
```

**性能注意**：线圈级函数是库中最耗时的调用（n_nodes=9 时 ~毫秒级，n_nodes=4 时 ~数十微秒）。在时间步进仿真中，通常每级每步调用一次，而非在紧密内循环中使用。

### use_cache 重载

每个丝级和线圈级函数都有一个接受末尾 `bool use_cache = false` 参数的重载：

```cpp
double mutual_inductance_filament(double radius_a, double radius_b,
                                  double separation, bool use_cache);
double mutual_inductance_gradient_filament(double radius_a, double radius_b,
                                           double separation, bool use_cache);
double mutual_inductance_coil(..., int n_nodes, bool use_cache);
double mutual_inductance_gradient_coil(..., int n_nodes, bool use_cache);
```

当 `use_cache` 为 `true` 时，函数会读写全局 4096 条 LRU 缓存。当相同的 `(radius_a, radius_b, separation)` 三元组被重复查询时，可以提升性能——典型场景是串行冷路径计算，如 `[M]` 矩阵初始化。

当 `use_cache` 为 `false`（默认值）时，缓存被完全绕过。在热路径（每步 `[M_I]` 更新）中使用此模式，因为缓存命中率低且需要线程安全。

不含 `use_cache` 的无参重载为了向后兼容而保留，其行为等同于 `use_cache = false`。

---

## LRU 缓存

```cpp
#include <coilgun/physics/cache.hpp>

template<typename Key, typename Value,
         std::size_t MaxSize = 4096,
         typename KeyHash = std::hash<Key>>
class coilgun::physics::LRUCache {
public:
    bool        get(const Key& key, Value& value);   // 命中返回 true
    void        put(const Key& key, const Value& value); // 插入或更新
    std::size_t size() const;                         // 当前条目数
    void        clear();                              // 清空所有条目
};
```

标准的 `unordered_map + list` LRU。非线程安全——每个调用方应持有自己的缓存实例。

由 `mutual_inductance_filament`（4096 条）和 `mutual_inductance_gradient_filament`（4096 条）内部使用。外部用户也可以自行实例化：

```cpp
LRUCache<std::string, Eigen::Vector3d, 1024> position_cache;
position_cache.put("home", {0.0, 0.0, 0.0});
Eigen::Vector3d pos;
if (position_cache.get("home", pos)) { /* 命中 */ }
```

---

## DrivingCoil（驱动线圈）

```cpp
#include <coilgun/components/driving_coil.hpp>

class coilgun::components::DrivingCoil {
public:
    DrivingCoil(double inner_radius, double outer_radius, double length,
                int turns, double resistivity, double wire_area,
                double fill_factor, double position = -1.0,
                bool force_exact_self_inductance = false);

    // 几何参数
    double inner_radius() const;     // ri (m)
    double outer_radius() const;     // re (m)
    double length() const;           // 轴向长度 (m)
    double mean_radius() const;      // (ri + re) / 2 (m)
    int    turns() const;            // N

    // 预计算电气属性
    double turns_density() const;    // nc = N / ((re - ri) × l)  (匝/m²)
    double resistance() const;       // 直流电阻 (Ω)，已考虑填充因子
    double self_inductance() const;  // 自感 (H)，通过 T(q,p) 查表或精确计算

    // 位置（仿真中用于平移）
    double position() const;         // 当前中心位置 (m)
    void   set_position(double x);   // 移动线圈中心至 x
};
```

所有电气属性在构造时预计算完成。`position` 字段用于在仿真中移动线圈而无需重新计算几何参数。

**构造参数说明**：

| 参数 | 说明 | 典型值 |
|-----------|-------------|---------------|
| `inner_radius` | 内绕组半径 | 0.01 m |
| `outer_radius` | 外绕组半径 | 0.03 m |
| `length` | 轴向绕组长度 | 0.05 m |
| `turns` | 总匝数 | 100–500 |
| `resistivity` | 参考温度下的导线电阻率 | `COPPER.resistivity_ref` |
| `wire_area` | 导体截面积 | 1e-6 m² (1 mm²) |
| `fill_factor` | 绕组填充因子 (0–1) | 0.5–0.8 |
| `position` | 初始中心位置 (m)，默认 = length/2 | 0.0 m |
| `force_exact_self_inductance` | 强制使用精确自感计算（跳过 T(q,p) 表）| false |

电阻公式已计入填充因子：仅 `fill_factor × 截面积` 部分参与导电。

---

## Armature（电枢）

```cpp
#include <coilgun/components/armature.hpp>

class coilgun::components::Armature {
public:
    Armature(double inner_radius, double outer_radius, double length,
             double resistivity, double material_density,
             double velocity, double mass,
             int m_axial, int n_radial, double position,
             physics::ArmatureMaterial material = physics::ArmatureMaterial::Aluminum,
             bool force_exact_self_inductance = false);

    // 几何参数
    double inner_radius() const;
    double outer_radius() const;
    double length() const;
    int    axial_filaments() const;     // m（轴向分段数）
    int    radial_filaments() const;    // n（径向分层数）
    int    total_filaments() const;     // m × n

    // 电流丝查询 — 1 基索引，与 NumericalModel 符号一致
    double filament_inner_radius(int j) const;         // j = 1..n
    double filament_outer_radius(int j) const;
    double filament_mean_radius(int j) const;          // 环中心线半径
    double filament_axial_position(int i) const;       // i = 1..m

    // 逐丝数据数组 — 长度为 m×n，行优先：(i=1,j=1), (i=1,j=2), ...
    const std::vector<double>& resistances() const;    // 每丝电阻 (Ω)
    const std::vector<double>& inductances() const;    // 每丝自感 (H)
    const std::vector<double>& masses() const;         // 每丝质量 (kg)

    // 运动状态
    double position() const;            // 当前中心位置 (m)
    double velocity() const;            // 当前速度 (m/s)
    double mass() const;                // 总质量（含载荷）(kg)
    void   update_position(double dx);  // 平移 dx 米
    void   set_velocity(double v);      // 设置新速度 (m/s)

    // 材料
    physics::ArmatureMaterial material() const;  // 热模式调度用
};
```

电枢是一个厚壁空心圆柱，离散为 m 个轴向段 × n 个径向层的电流丝。每根电流丝视为承载均匀电流的环；其电气属性（R、L）和质量在构造时预计算。

**构造参数说明**：

| 参数 | 说明 | 默认值 |
|-----------|-------------|---------|
| `inner_radius` | 内孔半径 (m) | — |
| `outer_radius` | 外半径 (m) | — |
| `length` | 轴向长度 (m) | — |
| `resistivity` | 材料电阻率 (Ω·m)，例如 `ALUMINUM.resistivity_ref` | — |
| `material_density` | 材料密度 (kg/m³) | — |
| `velocity` | 初始速度 (m/s) | — |
| `mass` | 含载荷总质量 (kg) — 各丝质量之和约为电枢材料质量；差值为载荷质量 | — |
| `m_axial` | 轴向分段数（典型 5–20） | — |
| `n_radial` | 径向分层数（典型 1–5） | — |
| `position` | 初始中心位置 (m) | — |
| `material` | 电枢材料：`Aluminum` 或 `Copper`（热模式下决定 cp(T) 和 beta）| `Aluminum` |
| `force_exact_self_inductance` | 强制使用精确纤维自感计算 | `false` |

**电流丝索引**：采用 1 基索引，与 NumericalModel.md 中的约定一致。`filament_axial_position(1)` 给出最左端轴向环的位置，`filament_axial_position(m)` 给出最右端。

**数据布局**：三个逐丝数组（`resistances()`、`inductances()`、`masses()`）为长度为 `m × n` 的扁平 vector，**行优先**排列。第 i 个轴向段、第 j 个径向层的索引为 `(i - 1) * n + (j - 1)`。

**示例**——一个 100g 铝质电枢加 20g 载荷，5 轴向 × 2 径向丝：

```cpp
Armature arm(
    0.005,                     // 内半径：5 mm 内孔
    0.025,                     // 外半径：25 mm
    0.08,                      // 长度：80 mm
    ALUMINUM.resistivity_ref,  // 2.82e-8 Ω·m
    ALUMINUM.density,           // 2700 kg/m³
    0.0,                       // 初始速度
    0.120,                     // 总质量 = 120g（100g 电枢 + 20g 载荷）
    5, 2,                      // m=5 轴向, n=2 径向 = 共 10 丝
    0.0                        // 初始位置
);

std::cout << "电流丝总数：" << arm.total_filaments() << "\n";               // 10
std::cout << "第(i=3,j=1)环半径：" << arm.filament_mean_radius(1) << " m\n";
std::cout << "第(i=3,j=1)环位置：" << arm.filament_axial_position(3) << " m\n";

// 通过扁平索引访问逐丝数据：
int idx = (3 - 1) * 2 + (1 - 1);   // (i=3, j=1) → 索引 4
std::cout << "电阻：" << arm.resistances()[idx] << " Ω\n";
std::cout << "质量：" << arm.masses()[idx] << " kg\n";
```

---

## 仿真器

`coilgun::simulation` 命名空间提供了一个完整的单级/多级线圈炮仿真器，支持可配置的激励源、时间积分策略和热耦合。

### 激励框架

```cpp
#include <coilgun/simulation/excitation.hpp>
```

| 类 | 描述 | 关键构造参数 |
|-------|-------------|---------------------|
| `Excitation` | 抽象基类 — `voltage()`、`advance(dt, I_coil)`、`finished()`、`reset()` | — |
| `CapacitorExcitation` | 电容放电（无续流二极管） | `(初始电压, 电容量)` |
| `CrowbarExcitation` | 电容放电（**含**续流二极管） | `(初始电压, 电容量)` |
| `WaveformExcitation` | 任意电压源 `V(t)` | `(V(t) 函数)` |

所有激励源均提供 `voltage()`、`advance(dt, I_coil)`、`finished()` 和 `reset()`。`CapacitorExcitation` 额外提供 `capacitance()`、`capacitor_voltage()` 和 `initial_voltage()`。`CrowbarExcitation` 通过 `diode_on()` 报告续流二极管状态。`WaveformExcitation` 支持通过 `set_end_time(t)` 设置可选的提前终止时间。

```cpp
auto cap  = std::make_unique<CapacitorExcitation>(450.0, 0.001);  // 450 V, 1000 μF
auto crow = std::make_unique<CrowbarExcitation>(450.0, 0.001);    // 含续流二极管
auto wfm  = std::make_unique<WaveformExcitation>(
    [](double t) { return 450.0 * std::exp(-t / 0.01); });
wfm->set_end_time(0.05);  // 在 50 ms 时停止波形

// 查询激励状态：
double v_cap  = crow->capacitor_voltage();  // 当前电容电压
double v_init = crow->initial_voltage();    // 初始电容电压 (450 V)
bool diode_on = crow->diode_on();           // 续流二极管是否导通？
crow->reset();                               // 恢复初始状态
```

### 时间步进器

```cpp
#include <coilgun/simulation/time_stepper.hpp>

// EulerStepper     — 前向 Euler（1 阶）
// RK4Stepper       — 经典 4 阶 Runge-Kutta
```

二者均为模板策略类。作为 `SingleStageSim` 或 `MultiStageSim` 的模板参数传入：

```cpp
SingleStageSim<EulerStepper> sim_euler(...);   // 较快，精度较低
SingleStageSim<RK4Stepper>   sim_rk4(...);     // 较慢，精度更高

MultiStageSim<EulerStepper> ms_euler(...);     // 多级同理
MultiStageSim<RK4Stepper>   ms_rk4(...);
```

### SimState

```cpp
#include <coilgun/simulation/single_stage_sim.hpp>

struct coilgun::simulation::SimState {
    Eigen::VectorXd currents;              // [N_fil+1], 索引 0 = 线圈电流, [1..N_fil] = 各丝
    double          arm_position = 0.0;    // m
    double          arm_velocity = 0.0;    // m/s
    Eigen::VectorXd filament_temperatures; // [N_fil], K（仅启用热耦合时有值）

    SimState& operator+=(const SimState& rhs);
    SimState& operator*=(double scalar);
};

SimState operator+(SimState lhs, const SimState& rhs);
SimState operator*(double scalar, SimState s);
```

`SingleStageSim::state()` 返回当前内部 `SimState` 的引用。算术运算符供内部步进算法使用。

### MultiStageState

```cpp
#include <coilgun/simulation/multi_stage_sim.hpp>

struct coilgun::simulation::MultiStageState {
    Eigen::VectorXd currents;              // [n_stages + N_fil], 先是各级线圈，后是各丝
    double          arm_position = 0.0;    // m
    double          arm_velocity = 0.0;    // m/s
    Eigen::VectorXd filament_temperatures; // [N_fil], K（仅启用热耦合时有值）

    MultiStageState& operator+=(const MultiStageState& rhs);
    MultiStageState& operator*=(double scalar);
};

MultiStageState operator+(MultiStageState lhs, const MultiStageState& rhs);
MultiStageState operator*(double scalar, MultiStageState s);
```

`MultiStageSim::state()` 返回当前内部 `MultiStageState` 的引用。

### 仿真结果

```cpp
#include <coilgun/simulation/sim_result.hpp>

struct SimStep {
    double time;                          // 时间
    double cap_voltage;                   // 电容电压
    double coil_current;                  // 线圈电流
    std::vector<double> filament_currents; // 各丝电流
    double arm_position;                  // 电枢位置
    double arm_velocity;                  // 电枢速度
    double force;                         // 有符号瞬时合轴向力 (N)
    std::vector<double> filament_temperatures; // 各丝温度（仅启用热耦合时）
};

struct SimSummary {
    double muzzle_velocity;    // 出口速度 (m/s)
    double total_time;         // 总耗时 (s)
    double max_force;          // 有符号瞬时力的峰值绝对值 (N)
    double peak_coil_current;  // 峰值线圈电流 (A)
    double efficiency;         // 效率 0–1
    int    step_count;         // 步数
};

struct SimResult {
    std::vector<SimStep> history;
    SimSummary           summary;
    SimResult            sampled(int every_n) const;  // 降采样，便于绘图
};

// 使用 sampled() 缩减导出数据：
auto sparse = result.sampled(100);  // 每 100 步保留 1 个
```

### 终止策略

```cpp
#include <coilgun/simulation/termination.hpp>

struct TerminationPolicy {
    int    max_steps             = 20000;
    int    velocity_decay_steps  = 5;      // 连续速度下降步数
    double accel_threshold       = 0.1;    // m/s² — 接近零加速度视为已停止
    bool   enable_velocity_check = true;
    bool   enable_bound_check    = false;
    double barrel_end_position   = 1.0;

    static TerminationPolicy defaults();
};

// 自定义策略：
TerminationPolicy pol = TerminationPolicy::defaults();
pol.max_steps = 50000;
sim.run(pol);
```

### SingleStageSim

```cpp
#include <coilgun/simulation/single_stage_sim.hpp>

template<typename StepperPolicy = EulerStepper>
class SingleStageSim {
public:
    SingleStageSim(
        components::DrivingCoil          coil,          // 拷贝传入
        components::Armature             armature,      // 拷贝传入
        std::unique_ptr<Excitation>      excitation,    // 移动传入
        double                           dt,
        bool                             enable_thermal = false
    );

    const SimStep&   step();                       // 前进一步
    const SimResult& run();                        // 运行至默认终止条件
    const SimResult& run(const TerminationPolicy& p); // 运行至自定义终止条件
    void             reset();                      // 重置为初始状态

    const SimResult& result()     const;   // 最近一次 run() 的结果
    const SimState&  state()      const;   // 内部状态（电流、位置、速度、温度）
    double           dt()         const;   // 时间步长
    int              step_count() const;   // 自构造或最后 reset() 以来的步数
};

// 以不同参数重新运行：
sim.run();
double v1 = sim.result().summary.muzzle_velocity;
sim.reset();          // 恢复初始状态（包括激励源）
sim.run();
double v2 = sim.result().summary.muzzle_velocity;  // 与 v1 完全相同
```

### OptimizationLevel（优化等级）

```cpp
#include <coilgun/simulation/multi_stage_sim.hpp>

enum class OptimizationLevel {
    Reference   = 0,   // 固定 9 点互感求积，无距离截断
    LookupTable = 1,   // 保留给组件级 T(q,p) 查表选择；当前运行时路径与 Reference 相同
    Full        = 2    // 距离截断 + 互感自适应 9/4 点求积
};
```

等级为累积叠加，每个等级包含所有较低等级的优化。

| 等级 | T(q,p) 查表 | 距离截断 | 自适应 GL 阶数 |
|:---:|:---:|:---:|:---:|
| 0 参考 | 否 | 否 | 否（固定 9 点） |
| 1 查表 | 是 | 否 | 否（固定 9 点） |
| 2 全优化 | 是 | 是 | 是（近=9，远=4） |

### 触发配置

```cpp
#include <coilgun/simulation/trigger_config.hpp>

enum class TriggerMode {
    Position,    // 当电枢中心越过指定位置时触发（m）
    TimeDelay    // 前一级触发后经过固定时间延迟触发（s）
};

struct TriggerConfig {
    TriggerMode mode;
    double      value;       // 触发位置 (m) 或时间延迟 (s)
};
```

第 0 级（stage 0）在 t=0 自动触发，无需 TriggerConfig。`trigger_configs[i-1]` 对应 stage `i`（i ≥ 1）。

`validate_trigger_config()` 对两种模式均拒绝无效模式和 NaN。`Position` 接受所有有限值，但拒绝负无穷；`TimeDelay` 接受有限非负值，并拒绝负值和负无穷。两种模式都允许正无穷，它是明确的终止性永不触发策略，表示不可到达的位置或延迟，而不是格式错误。有限的未触发 stage 即使在前一级完成后仍具备触发资格，因此会保持仿真运行，直到它触发或其他终止条件生效。只有每一级都已完成，或通过 `+infinity` 策略变为终止性不适用时，才会自动完成。触发位置在该级触发时的步前边界捕获。

### 多级仿真结果

```cpp
#include <coilgun/simulation/multi_stage_result.hpp>

struct StepSnapshot {
    double time;                    // 仿真时间 s
    double arm_position;            // 电枢中心位置 m
    double arm_velocity;            // 电枢速度 m/s
    double force;                   // 有符号瞬时合轴向力 N
    std::vector<double> filament_currents;       // 各纤维电流 A
    std::vector<double> filament_temperatures;   // 各纤维温度 K（仅热模式）
};

struct MultiStageStep {
    StepSnapshot          state;
    std::vector<double>   cap_voltages;     // [n_stages]，未触发级为 0
    std::vector<double>   coil_currents;    // [n_stages]，未触发级为 0
    std::vector<double>   stage_forces;     // [n_stages]，提交后的有符号瞬时各级力贡献 N
};

struct PerStageSummary {
    int    stage_index;           // 零基级号
    double trigger_time;          // 触发时刻 s
    double trigger_position;      // 触发步前边界的电枢位置 m
    double peak_current;          // 该级峰值电流 A
    double max_force;             // 该级有符号力贡献的峰值绝对值 N
    double energy_depleted;       // 电容消耗能量 J
    int    step_count_active;     // 该级活跃步数
};

struct MultiStageSummary {
    double                        muzzle_velocity;      // 出口速度 m/s
    double                        total_time;           // 仿真总时间 s
    double                        max_force;            // 有符号合力的峰值绝对值 N
    double                        peak_coil_current;    // 全局峰值电流 A
    double                        efficiency;           // 效率 0..1，E_kin / Σ(½C_i·U0_i²)
    int                           step_count;           // 总步数
    std::vector<PerStageSummary>  per_stage;            // 每级汇总
};

struct MultiStageResult {
    std::vector<MultiStageStep> history;   // 逐步状态记录
    MultiStageSummary           summary;   // 汇总统计
    MultiStageResult            sampled(int every_n) const;  // 降采样用于绘图
};
```

### MultiStageSim 类

```cpp
#include <coilgun/simulation/multi_stage_sim.hpp>

template<typename StepperPolicy = EulerStepper>
class MultiStageSim {
public:
    static constexpr int kMaxStages = 50;   // 最大级数上限

    MultiStageSim(
        std::vector<components::DrivingCoil>        coils,          // 拷贝传入
        components::Armature                         armature,      // 拷贝传入
        std::vector<std::unique_ptr<Excitation>>     excitations,   // 移动传入
        std::vector<TriggerConfig>                   trigger_configs,
        double                                       dt,
        bool                                         enable_thermal = false,
        OptimizationLevel                            opt_level = OptimizationLevel::Reference
    );

    const MultiStageStep&  step();     // 前进一步
    const MultiStageResult& run();     // 运行至默认终止条件
    const MultiStageResult& run(const TerminationPolicy& policy);
    void                    reset();    // 重置至初始状态

    const MultiStageResult& result()      const;
    const MultiStageState&  state()       const;
    double  dt()              const;   // 时间步长
    int     step_count()      const;   // 步数
    int     num_stages()      const;   // 已配置级数
};
```

**构造函数约束**：`coils.size() == excitations.size()`，`trigger_configs.size() == coils.size() - 1`，`coils.size() <= kMaxStages (50)`。

**系统说明**：ODE 维度为 `n_stages + 电枢纤维数`。未触发级（行列）设为恒等矩阵，保证矩阵非奇异。每级的 Crowbar 续流二极管自主管理，不需要额外矩阵操作。终止条件为所有级均 finished（diode ON 且电流衰减至零）。

**热模式**：当 `enable_thermal = true` 时，每个电枢纤维经历绝热焦耳加热（NumericalModel §6）。电阻每步根据温度相关电阻率更新。

**优化**：`opt_level` 参数控制运行时互感计算，独立于组件构造时确定的自感计算模式。`Reference` 和 `LookupTable` 当前使用相同运行时路径：固定 9 点 Gauss-Legendre 求积且无距离截断。`LookupTable` 仅为 API 兼容保留，不会回溯改变组件已经计算好的自感。`Full` 会跳过远处级，并对近处使用 9 点、远处使用 4 点求积。

---

## 完整示例

### 单级仿真

完整的单级仿真，含续流二极管、Euler 积分，无热耦合：

```cpp
#include <coilgun/coilgun.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::SingleStageSim;
    using coilgun::simulation::CrowbarExcitation;

    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    auto excitation = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(excitation), 1e-6, false);

    sim.run();

    const auto& r = sim.result();
    std::cout << "出口速度：" << r.summary.muzzle_velocity << " m/s\n";
    std::cout << "效率：    " << r.summary.efficiency * 100 << " %\n";
    std::cout << "步数：    " << r.summary.step_count << "\n";
    return 0;
}
```

`SingleStageSim` 自动处理完整的电路 ODE（线圈 + 电枢丝 + 互耦）、洛伦兹力积分、运动学更新和自动终止——全部通过所选的步进策略完成。

### 多级仿真

完整的两级仿真，含续流二极管和位置触发：

```cpp
#include <coilgun/coilgun.hpp>
#include <iostream>
#include <memory>
#include <vector>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::MultiStageSim;
    using coilgun::simulation::CrowbarExcitation;
    using coilgun::simulation::TriggerConfig;
    using coilgun::simulation::TriggerMode;

    // 构造两级线圈
    DrivingCoil coil1(0.01, 0.03, 0.05, 150,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.01, 0.03, 0.05, 150,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.10);

    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    std::vector<DrivingCoil> coils = {coil1, coil2};

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));

    // stage 1 在电枢到达 0.09m 时触发
    std::vector<TriggerConfig> triggers;
    triggers.push_back({TriggerMode::Position, 0.09});

    MultiStageSim<EulerStepper> sim(
        std::move(coils), arm, std::move(excs), triggers, 1e-6);

    sim.run();

    const auto& r = sim.result();
    std::cout << "出口速度：" << r.summary.muzzle_velocity << " m/s\n";
    std::cout << "效率：    " << r.summary.efficiency * 100 << " %\n";
    std::cout << "级数：    " << r.summary.per_stage.size() << "\n";
    for (const auto& ps : r.summary.per_stage) {
        std::cout << "  级 " << ps.stage_index
                  << ": 触发时刻=" << ps.trigger_time
                  << "s, 峰值电流=" << ps.peak_current << " A\n";
    }
    return 0;
}
```

`MultiStageSim` 自动处理完整的多级电路 ODE（级间互感耦合、反电动势、逐级自主 crowbar 二极管）、洛伦兹力积分、运动学更新、自动级触发和终止——全部通过所选的步进策略完成。

---

## 并行执行

仿真引擎使用 OpenMP 实现共享内存多核并行。每个时间步有三个计算循环被并行化：

- 线圈到电流丝互感计算（`M` 和 `dM/dx`）
- 洛伦兹力求和（归约）
- 电流丝温度更新（启用热模式时）

### 控制线程数

```cpp
#include <omp.h>
omp_set_num_threads(4);
```

或通过环境变量：

```sh
OMP_NUM_THREADS=4 ./your_simulation
```

默认线程数等于逻辑 CPU 数量。

### 线程安全

`mutual_inductance_filament` 和 `mutual_inductance_gradient_filament` 函数在 `use_cache = false`（默认值）时是线程安全的。全局 LRU 缓存仅在串行冷路径代码中访问（`[M]` 矩阵初始化、级间互感预计算）。

### CMake 集成

```cmake
find_package(OpenMP REQUIRED)
target_link_libraries(your_target PRIVATE coilgun OpenMP::OpenMP_CXX)
```

如果链接静态库 `libcoilgun.a`，还必须在你的目标中链接 OpenMP。

### 编译标志

所有 CMake 预设都包含 `-march=native`，用于 SIMD 指令集自动检测。Eigen 在编译时利用此标志为矩阵运算启用 AVX2/FMA/AVX512。

---

## GPU 加速 (CUDA)

本库提供可选的 GPU 加速后端 (`libcoilgun_cuda.a`)，将 4D Gauss-Legendre 互感积分卸载到 GPU。CPU 路径完全不变——GPU 类是独立的、直接替换的实现。

```sh
cmake --preset ninja-cuda-debug      # CUDA 配置
cmake --build --preset ninja-cuda-debug  # 构建
ctest --preset debug                 # 运行全部测试（CPU + GPU）
```

### 前置条件

| 需求 | 最低版本 | 备注 |
|------|----------|------|
| CUDA Toolkit | ≥ 12.8（Blackwell sm_120）；≥ 9.0（早期架构） | `nvcc` 必须在 PATH 中 |
| Boost.Math | ≥ 1.86 | 通过 CMake 自动获取 |
| NVIDIA GPU | 计算能力 ≥ 6.0（Pascal+） | FP64 原子操作支持所必需 |
| Host 编译器 | g++ ≥ 9 或同等 | 必须与 nvcc 检测到的一致 |

### 引入和链接

```cpp
#include <coilgun/coilgun_cuda.hpp>   // 所有 GPU 类 + CPU umbrella
```

```cmake
target_link_libraries(your_target PRIVATE coilgun_cuda coilgun CUDA::cudart)
```

`coilgun_cuda` 目标传递依赖 `coilgun`——两个库都需要。

---

### GpuBackend

```cpp
#include <coilgun/simulation/cuda/gpu_backend.hpp>

namespace coilgun::simulation::cuda {

struct GpuBackend {
    int     device_id         = 0;     ///< cudaSetDevice 目标。
    int     threads_per_block = 512;   ///< 4D 积分 kernel 的每 block 线程数。
    size_t  max_batch_sims    = 256;   ///< 批量仿真缓冲区的预分配上限。
    bool    enable_profiling  = false; ///< 保留 profiling 请求元数据；主机墙钟计时字段始终采集。不保证 NVTX。
    bool    use_persistent    = true;  ///< 遗留多级/批量包装器的持久化请求默认值。
    BackendMode backend        = BackendMode::Graph; ///< 未请求持久化时的后端。
};

}
```

| 字段 | 说明 | 默认值 |
|-------|---------|---------|
| `device_id` | 选择目标物理 GPU（多 GPU 系统相关）。 | `0` |
| `threads_per_block` | 4D GL 积分 kernel 请求的每 block 线程数。必须是不大于 512 的正数 2 的幂；1、128、256、512 均有效。 | `512` |
| `max_batch_sims` | `SimBatch` 中的最大仿真数。`SimBatch` 对超过该值的 `num_sims` 抛出 `std::invalid_argument`；该字段本身不负责分配缓冲区。 | `256` |
| `enable_profiling` | 为 true 时，在 `ExecutionReport::profiling_enabled` 中保留请求。主机墙钟计时类别独立于此标志始终采集。本构建不承诺 NVTX 标记，也不引入 NVTX 依赖。 | `false` |
| `use_persistent` | 遗留兼容请求。省略 `GpuMultiStageSim` 后端时使用 `multi_stage_default_backend()` 并选择 `Direct`；调用方显式提供 `use_persistent=true` 时仍有意请求 `Persistent`。在 `SimBatch` 中，只有当 `backend.backend == Auto` 时才读取该标志；显式后端模式对两种标志值保持相同的确定性和能力策略。 | `true` |
| `backend` | 后端元数据/配置。在 `GpuSingleStageSim` 中，`use_persistent=false` 时可选择 `Graph`、`Direct` 或 `Fallback`。在 `GpuMultiStageSim` 中，应使用构造函数的显式 `BackendMode` 参数，将有意的 `Graph`、`Direct` 或 `Fallback` 请求与遗留标志区分开。Graph 只捕获 mutual segment；CUDA 不可用或运行时失败时报告 `Fallback`。 | `Graph` |

故障注入控制不属于 `GpuExecutionConfig`。聚焦测试通过独立的内部 `detail::GpuEngineFaultInjection` 构造函数 seam 使用它们；普通执行配置只包含运行时策略和经过校验的 launch 设置。

**迁移状态**：`GpuEngine` 是当前执行核心。`GpuSingleStageSim`、`GpuMultiStageSim` 和 `SimBatch` 使用引擎契约进行后端选择、回滚、报告和部分 Graph 捕获。Graph 报告只表示 mutual segment 的图辅助执行，不表示完整多级流水线已捕获。

**示例**：

```cpp
coilgun::simulation::cuda::GpuBackend be;
be.device_id = 0;
be.threads_per_block = 256;   // 不大于 512 的任意正数 2 的幂
```

---

### GpuOptLevel

```cpp
#include <coilgun/simulation/cuda/gpu_backend.hpp>

namespace coilgun::simulation::cuda {

enum class GpuOptLevel {
    Standard   = 0,   ///< n_nodes=9，无距离截断。用于调试/验证。
    Full       = 1,   ///< 距离截断（>10× 线圈长度）+ 固定 n_nodes=9。生产环境。
    Aggressive = 2,   ///< FP32 integrand + FP64 reduction，距离截断，n_nodes=9。大规模参数扫描。
};

}
```

| 等级 | 距离截断 | GL 阶数 | 行为 |
|:-----:|:---------------:|:--------:|-----------|
| `Standard` | 否 | n_nodes = 9 | 所有 (stage, filament) 对都处理。kernel 使用 9 个 GL 节点（6561 个积分点）。最安全的设置。 |
| `Full` | 是 | n_nodes = 9 | 距离电枢超过 10× 线圈长度的 stage 被完全跳过。在典型距离下所有活跃 stage 都在范围内。 |
| `Aggressive` | 是 | n_nodes = 9 | 与 `Full` 相同但 integrand 使用 FP32 算术（含 FP32 AGM 椭圆积分）和 FP64 累加。消费级 GPU 上可获得 3-5× 加速。精度损失微小。 |

**注意**：GPU 后端不支持自适应 GL 阶数（n_nodes=4/9 混合）。GPU 上使用 n_nodes=4 会导致共享内存归约的非确定性浮点漂移（仅 8 个 warp 中的 3 个有实际工作）。CPU `OptimizationLevel::Full` 同时使用距离截断和自适应积分；GPU `GpuOptLevel::Full` 仅使用距离截断。

---

### 统一执行配置与规划

以下主机端类型会在创建任何 CUDA 资源前解析请求的执行契约：

```cpp
enum class BackendMode {
    Auto = 0, Graph = 1, Persistent = 2, Fallback = 3, Direct = 4
};
enum class SolverMode {
    Auto = 0, Eigen = 1, Batched = 2, CuSolver = Batched
};
enum class PrecisionMode { Standard = 0, Full = 1, Aggressive = 2 };
enum class ThermalMode { Auto = 0, Disabled = 1, Cpu = 2, Gpu = 3 };
enum class FallbackReason {
    None = 0, CapabilityUnavailable = 1, DeterminismRequired = 2,
    RuntimeFailure = 3, MetadataConflict = 4
};

struct GpuExecutionConfig {
    BackendMode backend = BackendMode::Auto;
    SolverMode solver = SolverMode::Auto;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Auto;
    bool enable_calibration = false;
    bool deterministic = false;
    int device_id = 0;
    int threads_per_block = 512;
    bool enable_profiling = false;
    void validate() const;
};

struct GpuCapability {
    bool supports_graph = true;
    bool supports_persistent = true;
    bool supports_batched_solver = true;
    bool supports_gpu_thermal = true;
    bool persistent_is_deterministic = false;
};

struct GpuExecutionPolicy {
    BackendMode requested_backend = BackendMode::Auto;
    SolverMode requested_solver = SolverMode::Auto;
    PrecisionMode requested_precision = PrecisionMode::Full;
    ThermalMode requested_thermal = ThermalMode::Auto;
    BackendMode backend = BackendMode::Fallback;
    SolverMode solver = SolverMode::Eigen;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Disabled;
    FallbackReason backend_fallback_reason = FallbackReason::None;
    FallbackReason solver_fallback_reason = FallbackReason::None;
    FallbackReason thermal_fallback_reason = FallbackReason::None;
};

class GpuExecutionPlanner {
public:
    static constexpr GpuExecutionPolicy plan(
        std::size_t n_stages, std::size_t n_filaments, std::size_t batch_size,
        bool thermal_enabled, const GpuCapability&, const GpuExecutionConfig&) noexcept;
};
```

约束和静态规则：

- `GpuExecutionPlanner::plan` 是纯主机代码。它使用维度、显式请求、确定性模式和提供的能力快照，不检查 CUDA 或计时。
- 显式 `Fallback` 是纯 CPU 契约，不应创建 CUDA context。主机规划器为了策略检查会保留独立请求的求解器/热模式；`GpuEngine` 在实际运行时将其归一化为 `Fallback + Eigen`，并在需要时归一化为 CPU 热模式。
- 当 `supports_graph == false` 时，显式 `Graph` 解析为带 `CapabilityUnavailable` 的 `Fallback`；运行时依据解析后的后端执行，也不会创建 context。当 graph 能力可用时，当前引擎可以在构造时升级请求，但只捕获/重放 mutual segment。矩阵组装、求解、力/状态编排和热更新仍在图外执行。
- 显式 `Persistent` 当前解析为安全回退，因为同步引擎没有常驻 kernel 所需的独立控制 stream。标记为非确定性的能力也不能满足确定性请求。
- 在编译 CUDA、目标设备可用且初始化成功时，`Direct` 使用直接 CUDA mutual pipeline。主机规划器对 `Auto` 保守地解析为 `Fallback`；`GpuEngine` 在运行时检测到设备后可以将 `Auto` 升级为 `Direct`。调用方若需要固定的 CUDA 后端，应显式请求 `Direct` 或 `Graph`。
- `SolverMode::Batched` 需要 `supports_batched_solver`；否则解析为带 `CapabilityUnavailable` 的 Eigen。`SolverMode::Auto` 仅在规划器的大工作负载规则下选择 Batched，否则选择 Eigen。
- `ThermalMode::Gpu` 需要 `supports_gpu_thermal`；否则解析为带 `CapabilityUnavailable` 的 CPU 热模式。`Auto` 仅在大工作负载且有能力支持时选择 GPU 热模式。
- 当前大工作负载规则为 `batch_size >= 8 || n_stages + n_filaments >= 128`。这是静态规则，不是性能保证。
- 当 `enable_calibration` 为 true 时，构造过程在工作区初始化后执行一次单位批次求解器校准。校准不推进物理状态，`reset()` 也不会重复校准；`ExecutionReport::calibrated` 是该一次性操作的元数据。
- 能力快照不等于设备可用性证明。CUDA 运行时枚举、设备选择、context 创建或分配仍可能失败，并会报告运行时回退。调用方必须检查解析后的报告。

`GpuCapability` 和规划器公开存在是为了支持确定性的主机测试和策略检查；它们不表示 CUDA 运行时检测已经完成。

枚举语义：

| 类型 | 取值 |
|---|---|
| `BackendMode` | `Auto` 延迟到运行时选择；`Graph` 请求 mutual segment 的图辅助捕获/重放；`Persistent` 请求常驻协议，但当前 `GpuEngine` 会回退；`Fallback` 为纯 CPU；`Direct` 直接启动 CUDA mutual segment。 |
| `SolverMode` | `Auto` 应用静态工作负载规则；`Eigen` 为主机 LDLT 路径；`Batched` 请求 CUDA 批量求解器。`CuSolver` 是 `Batched` 的别名。 |
| `PrecisionMode` | `Standard` 为无距离截断的 FP64，`Full` 为带生产截断的 FP64，`Aggressive` 使用 FP32 integrand 和 FP64 reduction。 |
| `ThermalMode` | `Auto` 应用工作负载/能力规则；`Disabled` 省略热更新；`Cpu` 和 `Gpu` 在支持时选择对应热路径。 |
| `FallbackReason` | `None` 表示无回退；`CapabilityUnavailable` 表示静态/运行时能力不支持；`DeterminismRequired` 表示拒绝非确定性选择；`RuntimeFailure` 表示 context/分配/捕获/步执行失败；`MetadataConflict` 表示保守解析或当前未实现的请求。 |

配置与策略字段：

| 类型 | 字段组 | 含义 |
|---|---|---|
| `GpuExecutionConfig` | `backend`、`solver`、`precision`、`thermal` | 请求的执行模式。 |
| `GpuExecutionConfig` | `enable_calibration`、`deterministic` | 请求构造期一次性校准和确定性策略约束。 |
| `GpuExecutionConfig` | `device_id`、`threads_per_block` | 运行时设备选择和经过校验的 launch 宽度。 |
| `GpuExecutionConfig` | `enable_profiling` | 复制到报告的元数据；计时独立采集，不表示启用 NVTX。 |
| `detail::GpuEngineFaultInjection` | `fail_after_mutual`、`fail_allocation`、`fail_device_initialization`、`fail_graph_capture` | 内部/测试专用构造函数 seam；不是普通执行配置。 |
| `GpuCapability` | `supports_graph`、`supports_persistent`、`supports_batched_solver`、`supports_gpu_thermal` | 仅供静态规划使用的调用方能力快照。 |
| `GpuCapability` | `persistent_is_deterministic` | Persistent 是否可满足 `deterministic=true`。 |
| `GpuExecutionPolicy` | `requested_*` | 为审计保留的原始请求模式。 |
| `GpuExecutionPolicy` | `backend`、`solver`、`precision`、`thermal` | 驱动资源创建和执行的解析后模式。 |
| `GpuExecutionPolicy` | `*_fallback_reason` | 各维度的静态解析原因。 |

---

### GpuSingleStageSim

```cpp
#include <coilgun/simulation/cuda/gpu_single_stage_sim.hpp>

namespace coilgun::simulation::cuda {

template<typename SP = EulerStepper>
class GpuSingleStageSim {
public:
    GpuSingleStageSim(
        components::DrivingCoil          coil,          // 拷贝
        components::Armature             armature,      // 拷贝
        std::unique_ptr<Excitation>      excitation,    // 移动传入
        double                           dt,
         bool                             enable_thermal = false,
         GpuOptLevel                      opt_level = GpuOptLevel::Full,
         const GpuBackend&                backend = {},
         BackendMode                      explicit_backend = BackendMode::Auto
    );

    const SimStep&   step();
    const SimResult& run();
    const SimResult& run(const TerminationPolicy& p);
    void             reset();

    const SimResult& result()     const;
    const SimState&  state()      const;
    double  dt()         const;
    int     step_count() const;
    const ExecutionReport& execution_report() const;
    std::vector<double> filament_resistances() const;
};

} // namespace
```

**构造参数**：

| 参数 | 说明 | 典型值 |
|-----------|-------------|---------------|
| `coil` | 驱动线圈几何（拷贝） | — |
| `armature` | 电枢几何与丝元离散化（拷贝） | — |
| `excitation` | 激励源（移动传入） | `CrowbarExcitation(450.0, 0.001)` |
| `dt` | 固定时间步长（s） | `1e-6` |
| `enable_thermal` | 启用绝热丝元升温（CPU 侧） | `false` |
| `opt_level` | GPU 优化等级 | `GpuOptLevel::Full` |
| `backend` | GPU 后端配置 | `{}`（默认值） |
| `explicit_backend` | 可选的构造函数级后端覆盖。`Auto` 保留 `backend`；其他值覆盖 `backend.backend` 和 `use_persistent`。 | `BackendMode::Auto` |

当激励源为空、`dt` 非有限或不为正、`device_id` 为负、`threads_per_block` 无效、`max_batch_sims == 0`、激励电压非有限、几何/状态维度无效或 stage 电压非有限时，构造函数抛出 `std::invalid_argument`。启用 CUDA 时还会验证目标设备存在，并在分配资源前确认该设备仍被选中。`GpuSingleStageSim<RK4Stepper>` 为保持源码兼容仍可构造，但由于 RK4 不受支持，`step()` 抛出 `std::logic_error`。

**方法**：

| 方法 | 行为 |
|--------|-----------|
| `step()` | 推进一个 Euler 时间步。`Direct` 直接执行现有 CUDA mutual segment；`Graph` 捕获/重放该 mutual segment，然后在图外执行现有同步矩阵、Eigen LDLT、力/状态和热阶段。`Fallback` 在 CPU 上执行完整时间步。在 `Full`/`Aggressive` 模式下，当 `abs(电枢位置 - coil.position()) > 10 * coil.length()` 时，stage-电枢互感和力被截断，但驱动电路仍保持活动。返回新记录的 `SimStep`。 |
| `run()` | 运行至默认终止条件。阻塞直到激励结束（续流二极管导通且电流衰减完毕）。 |
| `run(policy)` | 运行至自定义 `TerminationPolicy`。 |
| `reset()` | 恢复至初始状态。所有电流归零，电枢位置/速度恢复至构造时值，重置激励，清空结果历史。 |
| `result()` | 最后一次 `run()` 之后的结果。包含 `history`（SimStep 向量）和 `summary`（SimSummary）。 |
| `state()` | 当前内部状态——实时引用。 |
| `dt()` | 固定时间步长。 |
| `step_count()` | 自构造或上次 `reset()` 以来的步数。 |
| `execution_report()` | 已解析的后端/求解器/热模型策略及累计诊断。回退时 `gpu_executed` 为 false；只有完整 CUDA 物理步成功提交后才为 true。GPU/求解器/热/传输计时、图重建、回退次数、校准、条件数估计和回退原因在 reset 后仍为累计值。 |
| `filament_resistances()` | 返回当前丝元电阻向量的值拷贝。启用热模型时，构造过程从电枢参考电阻初始化该向量；焦耳热会更新它，`reset()` 会恢复参考值。拷贝在包装器后续修改后仍有效，不暴露设备指针。 |

`SimStep::force` 是使用提交后的电流和在步前位置边界计算并缓存的互感梯度记录的有符号力。位置更新不会再次计算梯度；这是单级、多级和批量包装器共同采用的“步后电流/步前位置”约定。Euler 速度更新使用从步前电流计算的力。

`ExecutionReport` 字段语义：

| 字段 | 含义 |
|---|---|
| `requested_backend`、`requested_solver`、`requested_precision`、`requested_thermal` | 为审计保留的原始请求。 |
| `backend`、`solver`、`precision`、`thermal` | 当前解析后的执行模式。 |
| `gpu_executed` | 至少一个完整 CUDA 后端物理步成功提交的累计证明；不是能力/请求标志。 |
| `calibrated`、`precision_fallback`、`metadata_conflict` | 一次性校准完成状态及累计策略/报告诊断。 |
| `graph_rebuild_count`、`fallback_count` | `graph_rebuild_count` 只统计成功捕获的新 CUDA Graph 变体；主机变体选择、Direct/Fallback 构造、缓存命中和失败捕获都不增加该值。`fallback_count` 统计回退事件。 |
| `gpu_time_ms` | 成功 CUDA 后端物理流水线的累计主机墙钟时间，包括传输、CPU 矩阵组装/Eigen 求解和主机编排；绝不是仅设备 kernel 时间。 |
| `solver_time_ms`、`thermal_time_ms` | 实际 CPU 或 CUDA 后端路径中对应阶段的累计主机墙钟时间。 |
| `transfer_time_ms` | 同步主机/设备拷贝消耗的累计主机墙钟时间。 |
| `max_condition_estimate` | 已记录的最大求解器条件数估计/校准诊断。 |
| `fallback_reason` | 保留的最新人类可读回退消息。 |
| `static_fallback_reason`、`runtime_fallback_reason` | 机器可读的规划期和运行时原因。 |
| `device_id`、`threads_per_block` | 从配置复制的已校验执行设置。 |
| `profiling_enabled` | 仅为 profiling 请求元数据。计时独立采集，不表示存在 NVTX。 |

**所有权与回退**：包装器通过 RAII 持有一个 `GpuEngine`。`Excitation` 仍是移动传入的 `std::unique_ptr`；公开包装器不持有裸设备指针。`GpuEngine` 负责 CUDA context、求解器、图、热工作区和设备缓冲区的生命周期，并在初始化失败时释放已经部分分配的资源。CUDA 运行时枚举或设备选择失败、CUDA 不可用、context/分配失败或运行时流水线失败时，引擎在需要时恢复完整步前状态，初始化 Eigen 求解器，在 CPU 上执行完整时间步，并锁定后续步骤继续 CPU 回退。枚举/选择/context/分配/流水线错误会保留非空原因并设置 `runtime_fallback_reason=RuntimeFailure`；无设备/驱动不足使用 `CapabilityUnavailable`。成功 CUDA 步和 CPU 回退使用同一 Euler 物理流水线；但互感计算中的 CPU/GPU 浮点差异仍可能影响长时间运行结果。

**积分契约**：支持 `GpuSingleStageSim<EulerStepper>`。本次迁移中明确不支持 `GpuSingleStageSim<RK4Stepper>`：`step()` 会抛出 `std::logic_error`，不会静默改用 Euler。真正的四阶段 RK4 GPU 对齐被推迟，本 API 不宣称 RK4 一致性。

**重置与诊断**：`reset()` 清除仿真状态、激励状态、已完成步数和结果历史，然后在步边界重新选择图变体。`ExecutionReport` 是执行审计信息：回退次数、计时、图重建次数、校准状态、最大条件数估计、回退原因和 `gpu_executed` 均为累计字段，reset 后保留。启用校准时，构造过程会在工作区初始化后执行一次单位批次求解器校准；校准不会推进物理状态，reset 也不会重复校准。

**从 HEAD 迁移**：

HEAD 中的旧声明和默认值：

```cpp
struct GpuBackend {
    int device_id = 0;
    int threads_per_block = 512;
    size_t max_batch_sims = 256;
    bool enable_profiling = false; // 声称 NVTX range 标记
    bool use_persistent = true;
};
```

当前声明和默认值：

```cpp
struct GpuBackend {
    int device_id = 0;
    int threads_per_block = 512;
    size_t max_batch_sims = 256;
    bool enable_profiling = false; // 元数据和主机计时；不保证 NVTX
    bool use_persistent = true;
    BackendMode backend = BackendMode::Graph;
};
```

| 契约 | HEAD | 当前 E2 契约 |
|---|---|---|
| `GpuBackend` 声明/默认值 | `device_id=0`、`threads_per_block=512`、`max_batch_sims=256`、`enable_profiling=false`、`use_persistent=true`；没有 `backend` 字段 | 共享 `GpuBackend` 保留 `use_persistent=true`，并新增 `backend=BackendMode::Graph`。省略 `GpuMultiStageSim` 后端时使用 `multi_stage_default_backend()` 并请求 `Direct`；显式提供 `use_persistent=true` 仍表示有意请求 `Persistent`。显式构造模式覆盖该映射，可请求 `Graph`、`Direct`、`Fallback` 或 `Persistent`。`enable_profiling` 是元数据/计时请求，不再声称 NVTX。 |
| profiling/计时 | `enable_profiling` 声称 NVTX range，`gpu_time_ms` 容易被理解为设备时间 | `profiling_enabled` 只记录元数据。`gpu_time_ms` 是成功 CUDA 后端物理流水线的主机墙钟时间，包括传输和主机编排，不是仅设备时间。 |
| RK4 | 包装器执行所选的 `StepperPolicy`，包括遗留的 RK4 模板实例 | `GpuSingleStageSim<RK4Stepper>` 仍可构造，但 `step()` 抛出 `std::logic_error`，不会静默执行 Euler。 |
| 报告访问 | 没有统一包装器报告访问器 | `execution_report()` 返回累计 `ExecutionReport`，包括请求/解析模式、回退原因、计时、校准元数据和 `gpu_executed`。 |
| 热电阻访问 | 包装器没有更新后丝元电阻的访问器 | `filament_resistances()` 返回物理电阻向量的值拷贝；构造/重置时使用电枢参考值，焦耳热后更新；不暴露裸设备指针。 |
| 校验 | 遗留构造可能延后或遗漏边界校验 | 无效公开配置、空激励源、无效维度/几何/状态、非有限电压、无效 launch 设置以及负数/越界的配置设备 ID 抛出 `std::invalid_argument`。真实 CUDA 枚举/设备选择失败、无设备/驱动不足检测以及 context/分配/流水线失败则使用诚实且锁定的 CPU 回退；运行时错误设置 `runtime_fallback_reason=RuntimeFailure`，运行时能力不可用则使用 `CapabilityUnavailable`。 |
| 所有权 | 遗留包装器直接持有 adaptor/持久化资源 | 包装器通过 `std::unique_ptr` 持有 `GpuEngine`；`Excitation` 仍移动传入；公共包装器不暴露裸设备指针。引擎使用 RAII，包括部分分配清理。 |
| 回退 | 持久化和设备失败可能从请求或遗留路径推断 | 应检查解析后的 `backend`、`solver`、`gpu_executed`、`static_fallback_reason`、`runtime_fallback_reason` 和 `fallback_reason`。显式 Fallback 和能力不支持的 Graph 都是 CPU-only，且不创建 CUDA context。运行时可用性仍可能在能力快照后强制回退；配置了无效设备时则会在分配前拒绝。 |
| CPU/GPU 对齐 | 遗留 GPU 行为没有暴露迁移后的引擎契约 | 成功 GPU 路径和 CPU 回退共享 Euler 物理流水线，并由聚焦测试比较，但 CUDA/CPU 浮点差异仍可能存在。Graph 当前只捕获 mutual；矩阵、Eigen 求解器、力/状态和热阶段仍在图外。 |

迁移示例：

```cpp
// HEAD：persistent 是隐式默认值。
GpuBackend old_style{};

// E1：Graph 是显式默认值；构造/步进后检查报告。
GpuBackend graph;
graph.backend = BackendMode::Graph;
graph.use_persistent = false;
GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(excitation), 1e-6,
                                   false, GpuOptLevel::Full, graph);
sim.step();
const auto& report = sim.execution_report();
if (report.gpu_executed) {
    // report.gpu_time_ms 是已提交 CUDA 流水线的主机墙钟时间。
}

// 显式 CPU-only 迁移。
GpuBackend cpu;
cpu.backend = BackendMode::Fallback;
GpuSingleStageSim<EulerStepper> cpu_sim(coil, arm, std::move(cpu_excitation), 1e-6,
                                        false, GpuOptLevel::Full, cpu);
```

**完整示例**：

```cpp
#include <coilgun/coilgun_cuda.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::GpuSingleStageSim;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;

    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6);

    sim.run();
    std::cout << "出口速度: " << sim.result().summary.muzzle_velocity << " m/s\n";
    std::cout << "效率:     " << sim.result().summary.efficiency * 100 << " %\n";
    return 0;
}
```

构造和步进接口有意保持接近 `SingleStageSim`，但 GPU 包装器额外提供 `execution_report()`，并且当前只支持 Euler。CPU RK4 仿真在真正 GPU RK4 加入前必须继续使用 CPU 实现；GPU 包装器不宣称 RK4 对齐。

---

### GpuMultiStageSim

```cpp
#include <coilgun/simulation/cuda/gpu_multi_stage_sim.hpp>

namespace coilgun::simulation::cuda {

template<typename SP = EulerStepper>
class GpuMultiStageSim {
public:
    static constexpr int kMaxStages = 50;

    GpuMultiStageSim(
        std::vector<components::DrivingCoil>     coils,          // 拷贝
        components::Armature                      armature,      // 拷贝
        std::vector<std::unique_ptr<Excitation>>  excitations,   // 移动传入
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal = false,
        GpuOptLevel                               opt_level = GpuOptLevel::Full,
         const GpuBackend&                         backend = multi_stage_default_backend(),
        BackendMode                               explicit_backend = BackendMode::Auto
    );

    ~GpuMultiStageSim();

    GpuMultiStageSim(const GpuMultiStageSim&) = delete;
    GpuMultiStageSim& operator=(const GpuMultiStageSim&) = delete;

    const MultiStageStep&  step();
    const MultiStageResult& run();
    const MultiStageResult& run(const TerminationPolicy& p);
    void                    reset();

    const MultiStageResult& result()      const;
    const MultiStageState&  state()       const;
    double  dt()              const;
    int     step_count()      const;
    int     num_stages()      const;
    const ExecutionReport& execution_report() const;
    bool                    graph_assisted() const;
    std::vector<double> filament_resistances() const;
    std::vector<double> mutual_inductances() const;
    std::vector<double> mutual_gradients() const;
};

} // namespace
```

**构造约束**：`coils.size() == excitations.size()`，`trigger_configs.size() == coils.size() - 1`，`coils.size() ≤ kMaxStages (50)`，`dt` 为有限正数，每个激励源非空且电压有限，每个触发配置均通过 `validate_trigger_config()`，后端设置有效。违反任一则抛出 `std::invalid_argument`，包括无效的显式后端枚举值。有限位置和有限非负延迟始终保持触发资格；`+infinity` 是唯一有效的终止性永不触发值。

**构造参数**：

| 参数 | 说明 | 典型值 |
|-----------|-------------|---------------|
| `coils` | 每级驱动线圈几何（拷贝） | — |
| `armature` | 电枢几何与丝元离散化（拷贝） | — |
| `excitations` | 每级一个 `Excitation`（移动传入） | 各 `CrowbarExcitation(450.0, 0.001)` |
| `trigger_configs` | 除第 0 级外每级一个 `TriggerConfig` | `{Position, 0.09}` |
| `dt` | 固定时间步长（s） | `1e-6` |
| `enable_thermal` | 启用绝热丝元升温（CPU 侧） | `false` |
| `opt_level` | GPU 优化等级 | `GpuOptLevel::Full` |
| `backend` | GPU 后端配置。省略时使用 `multi_stage_default_backend()`（`use_persistent=false`、`backend=Direct`）。显式提供 `use_persistent=true` 仍表示有意请求 Persistent。 | `multi_stage_default_backend()` |
| `explicit_backend` | 可选的显式后端覆盖。`Auto` 保留所提供的遗留 `use_persistent` 行为；其他值覆盖 `use_persistent` 和 `backend.backend`。 | `BackendMode::Auto` |

**方法**——签名与 CPU 版 `MultiStageSim` 完全一致：

| 方法 | 行为 |
|--------|-----------|
| `step()` | 推进一个 Euler 时间步。已触发且未完成的 stage 始终参与线圈电路求解。在 `Full`/`Aggressive` 模式下，距离截断只作用于 stage-电枢互感项和力；stage 回到范围内前这些项为零。`Direct` 发射 CUDA mutual segment；`Graph` 仅捕获/重放该 mutual segment，并在图外执行矩阵组装、Eigen 求解、力/状态更新和热更新；`Fallback` 在 CPU 上执行完整步。`GpuMultiStageSim<RK4Stepper>::step()` 抛出 `std::logic_error`，不会静默执行 Euler。 |
| `run()` / `run(policy)` | 运行至终止条件。只有当每一级都已完成或通过 `+infinity` 变为终止性不适用时才会自动完成；有限的延时/位置触发级会保持运行，直到触发。显式最大步数、速度衰减和可选位置边界仍可终止具有可触发 stage 的运行。速度衰减判定使用最新提交的步后记录力计算加速度，与 CPU `compute_force(state_)` 一致，而不是积分时施加的步前力。 |
| `reset()` | 恢复所有状态至初始条件。 |
| `result()` / `state()` / `dt()` / `step_count()` / `num_stages()` | 查询。 |
| `execution_report()` | 读取请求/解析后的后端、`gpu_executed`、部分 Graph 计数、回退原因、计时和累计诊断。`gpu_executed` 是完整 CUDA 物理步已提交的证明，不只是请求或能力标志。 |
| `graph_assisted()` | 只有在完整 CUDA 后端步解析为 `BackendMode::Graph` 后才为 true；它表示 mutual segment 的图辅助，不表示完整流水线捕获。 |
| `filament_resistances()` | 返回当前丝元电阻向量的值拷贝。热模式构造时使用电枢参考值，热更新会改变它，`reset()` 恢复参考值。拷贝独立于包装器后续修改。 |
| `mutual_inductances()` / `mutual_gradients()` | 返回当前 stage/丝元 M1 和 dM1 缓冲区的值拷贝，使用引擎坐标。用于数值对齐检查，不暴露原始设备指针。 |

`MultiStageStep::state.force` 是状态更新后记录的有符号瞬时合轴向力。`stage_forces` 是对应的有符号瞬时各级力贡献，使用提交后的电流和在步前位置边界计算并缓存的互感梯度。如果激励在该步推进时完成某一级，则该级提交的力贡献为零，与 CPU 参考实现一致。用于推进速度的是单独保留的步前施加力，不是公开记录力；速度衰减终止判定会明确使用已提交的记录力。`PerStageSummary::max_force` 是该 stage 自身有符号记录力贡献的峰值绝对值；`MultiStageSummary::max_force` 是有符号记录合力历史的峰值绝对值。触发位置在实际触发的步前边界捕获，不再从后续带时间戳的历史样本反推。

**内部实现**：通过 `GpuEngine` 计算所有活跃 stage-丝元对的 M1/dM1 矩阵。Graph 模式只捕获/重放该 mutual-inductance segment；丝元/级间矩阵组装、Eigen LDLT 求解、力/状态编排、位置/速度更新、激励推进、触发处理和热更新均在图外执行。ODE 维度为 `n_stages + N_filaments`；未激活 stage 使用 mask/单位矩阵语义。必须检查 `execution_report().backend` 和 `gpu_executed`，以区分 Graph 辅助 CUDA 执行与 CPU 回退。

**E2 迁移契约**：

| 关注点 | 遗留/当前行为 | E2 行为 |
|---|---|---|
| `use_persistent` | `GpuBackend` 默认值仍为 `true` 以保持源码兼容；旧多级调用方使用 `false` 选择逐对 direct/fallback 路径。 | 省略 `GpuMultiStageSim` 后端时使用 `multi_stage_default_backend()` 并请求 `Direct`。显式提供 `use_persistent=false` 时请求 `Direct`（除非 `backend=Fallback`）；显式提供 `use_persistent=true` 时请求 `Persistent`，并可诚实回退。提供 `explicit_backend=Graph`、`Direct`、`Fallback` 或 `Persistent` 时显式覆盖且在 `requested_backend` 中可见。 |
| Graph | Graph 请求可能被理解为完整流水线。 | Graph 明确是图辅助：CUDA Graph 只捕获 mutual-inductance。`graph_assisted()` 和 `execution_report()` 标识实际成功的部分 Graph 执行；CPU 矩阵/求解/力/状态/热处理仍在图外。解析为 `Fallback` 且 `gpu_executed=false` 不算 GPU 执行。 |
| RK4 | CPU `MultiStageSim<RK4Stepper>` 执行真正的四阶段 RK4。 | `GpuMultiStageSim<RK4Stepper>` 为保持源码兼容仍可构造，但 `step()` 抛出 `std::logic_error`。不会实现伪 RK4，也不会静默替换为 Euler。存在真正的分阶段引擎 API 前，应迁移到 CPU `MultiStageSim<RK4Stepper>`。 |
| 生命周期 | 遗留 adaptor/持久化资源由包装器专门管理。 | `GpuEngine` 通过 RAII 管理 CUDA context、图/cache、求解器、热缓冲区和状态。`execution_report()` 仍是包装器持有的引用；`filament_resistances()` 返回独立值拷贝。reset 恢复仿真状态但保留累计报告诊断。 |
| 异常 | 一些无效输入可能延迟到下层才发现。 | 构造函数对级数/激励/触发数量不一致、空或过大的级数、非正/非有限 `dt`、空/非有限激励、无效后端模式或 launch 设置、无效几何/状态维度抛出 `std::invalid_argument`。CUDA 运行时不可用/失败则报告锁定的 CPU 回退，不与参数校验混淆。 |
| 对齐 | GPU-vs-CPU 测试可能在两条路径都是 CPU 回退时通过。 | CUDA 运行时可用时，主要非退化对比必须要求 `gpu_executed=true` 和预期解析后端，并比较电流、位置、速度、力、stage 输出以及热/电阻状态。另有独立显式 fallback 测试断言 CPU-only 执行。 |

**E2 迁移示例**：

```cpp
// 遗留多级语义：false 选择 direct/逐对路径。
GpuBackend legacy_direct;
legacy_direct.use_persistent = false;
GpuMultiStageSim<EulerStepper> direct(
    coils, armature, std::move(excitations), triggers, 1e-6,
    false, GpuOptLevel::Standard, legacy_direct); // requested Direct

// 有意请求 Graph：必须使用显式模式与遗留 false 标志区分。
// 它仅是图辅助执行。
GpuBackend graph_backend;
graph_backend.use_persistent = false;
GpuMultiStageSim<EulerStepper> graph(
    coils, armature, std::move(graph_excitations), triggers, 1e-6,
    false, GpuOptLevel::Standard, graph_backend, BackendMode::Graph);
graph.step();
const auto& report = graph.execution_report();
if (report.gpu_executed && graph.graph_assisted()) {
    // 只有 mutual-inductance segment 被捕获/重放。
}

// RK4 迁移仍使用 CPU。
MultiStageSim<RK4Stepper> cpu_rk4(
    coils, armature, std::move(rk4_excitations), triggers, 1e-6);
```

**完整示例**：

```cpp
#include <coilgun/coilgun_cuda.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::GpuMultiStageSim;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;
    using coilgun::simulation::TriggerConfig;
    using coilgun::simulation::TriggerMode;

    DrivingCoil c1(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    std::vector<DrivingCoil> coils = {c1, c2};
    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};

    GpuMultiStageSim<EulerStepper> sim(
        std::move(coils), arm, std::move(excs), triggers, 1e-6);

    sim.run();
    auto& r = sim.result();
    std::cout << "出口速度: " << r.summary.muzzle_velocity << " m/s\n";
    std::cout << "效率:     " << r.summary.efficiency * 100 << " %\n";
    for (auto& ps : r.summary.per_stage)
        std::cout << "  第" << ps.stage_index << "级"
                  << ": I_peak=" << ps.peak_current << " A\n";
    return 0;
}
```

**与 CPU `MultiStageSim` 的区别**：

| 方面 | CPU | GPU |
|--------|-----|-----|
| M/dM 计算 | OpenMP 并行化 | CUDA kernel（每对一个 block） |
| 优化等级 | `OptimizationLevel`（3 级） | `GpuOptLevel`（3 级） |
| 自适应 GL 阶数 | n_nodes=9/4（距离决定） | 仅 n_nodes=9 |
| 温度更新 | OpenMP 并行化 | CPU 串行循环 |
| 力求和 | OpenMP 归约 | 串行循环 |
| T(q,p) 查表 | 表范围内使用 | 仅构造时（非 GPU） |

---

### SimBatch

```cpp
#include <coilgun/simulation/cuda/sim_batch.hpp>

namespace coilgun::simulation::cuda {

template<typename SP = EulerStepper>
class SimBatch {
public:
    static constexpr int kMaxStages = 50;

    SimBatch(
        std::vector<components::DrivingCoil> coils,
        components::Armature                  armature,
         int                                   num_sims,
         double                                dt,
         const GpuBackend&                     backend = {},
         BackendMode                           explicit_backend = BackendMode::Auto);

    ~SimBatch();

    SimBatch(const SimBatch&) = delete;
    SimBatch& operator=(const SimBatch&) = delete;

    void set_excitations(
        int sim_id,
        std::vector<std::unique_ptr<Excitation>> excitations,
        std::vector<TriggerConfig>               trigger_configs);

    void run();
    void run(const TerminationPolicy& policy);

    const MultiStageResult& result(int sim_id) const;
    int num_sims() const;
    const ExecutionReport& execution_report() const;
    bool graph_assisted() const;
};

} // namespace
```

`SimBatch` 是**参数扫描**的容器——运行共享相同线圈和电枢几何但激励参数（电压、电容）和/或触发位置不同的多组仿真。

**执行契约**：`SimBatch` 是由 `GpuEngine` 驱动的包装器。一个引擎以 `B = num_sims` 构造；所有物理缓冲区都使用固定的行主序布局，包括 `[B][S][F]` 的互感/梯度数组和 `[B][S+F]` 的电流行。active mask 会原地冻结已完成的行。行不会压缩或重新编号，因此较短的仿真先结束时，`result(sim_id)` 仍保持稳定索引。这个固定容量策略是有意设计的：步骤之间不会压缩活跃行。活跃列表压缩阈值在 E3 中明确推迟且尚未实现；物理索引和结果索引保持稳定。

所有仿真必须共享**完全相同**的线圈几何和丝元离散化（`m × n`）。每个仿真的激励源、触发配置和 stage 电压通过 `set_excitations()` 提供。电路 mask 选择 stage 是否参与；互感 mask 还选择 stage-电枢互感和力。无论最终解析为何种后端（包括 `Fallback`），远距离 stage 都使用统一截断规则 `abs(armature_position - coil.position()) <= 10 * coil.length()`；超出范围时其互感和记录力项为零。

`SimBatch` 仅支持 `EulerStepper`。`SimBatch<RK4Stepper>` 会被明确拒绝：构造和 `run()` 抛出 `std::logic_error`，不宣称 RK4 对齐，也不会静默替换为 Euler。

**构造参数**：

| 参数 | 说明 |
|-----------|-------------|
| `coils` | 共享驱动线圈几何（拷贝）。所有仿真相同。 |
| `armature` | 共享电枢几何和丝元离散化（拷贝）。 |
| `num_sims` | 仿真数量。必须为正数且 ≤ `GpuBackend` 中的 `max_batch_sims`。 |
| `dt` | 固定时间步长（s），所有仿真共享。 |
 | `backend` | GPU 后端配置。构造函数的显式 `explicit_backend` 优先级最高；否则显式 `backend.backend` 对 `Graph`、`Direct`、`Fallback` 和 `Persistent` 生效；只有两者都为 `Auto` 时才读取遗留标志。解析后的模式和回退原因可从 `execution_report()` 获取。 |
 | `explicit_backend` | 可选的构造函数级后端覆盖。`Auto` 保留 `GpuBackend` 的解析结果；其他值覆盖 `backend.backend` 和 `use_persistent`。 |

**方法**：

| 方法 | 行为 |
|--------|-----------|
| `set_excitations(sim_id, excitations, triggers)` | 为仿真 `sim_id` 配置激励源和触发设置。必须在 `run()` 之前对每个仿真调用。 |
| `run()` / `run(policy)` | 同时运行所有仿真。每个仿真独立运行——最先到达终止条件的仿真不阻塞其他仿真。首次调用前必须配置每一行：`SimBatch` 的一次运行是单次使用的，后续调用不会重置或继续已完成的历史。 |
| `result(sim_id)` | 获取特定仿真的 `MultiStageResult`。 |
| `num_sims()` | 本批次中的仿真数。 |
| `execution_report()` | 返回共享引擎请求/解析后的后端和求解器、回退诊断、计时元数据以及 `gpu_executed`。 |
| `graph_assisted()` | 只有在完整 CUDA 步解析为 `Graph` 且实际执行 CUDA 后才为 true。Graph 辅助仅覆盖互感 segment。 |

**求解器和回退**：包装器请求 `SolverMode::Auto`，允许 `GpuExecutionPlanner` 在批量较大或维度较大时选择 `Batched`。如果 CUDA context 或批量能力不可用，引擎报告 `Eigen` 并执行 CPU 回退。`SimBatch` 不会因为请求了 `Batched` 就声称 cuSOLVER 已完成；应检查 `execution_report().solver`、`backend`、`gpu_executed` 和回退字段。

**步进和历史语义**：每一步引擎执行前，`SimBatch` 处理异构触发器、更新电路/互感 mask，并保存步前电流以及步前位置边界的互感梯度缓存。引擎使用这些值进行物理 Euler 更新。步进后，各激励源推进，已完成 stage 被 mask；记录的每-stage 力使用 post-step 电流和该步前位置边界的梯度缓存重新计算。这与 `GpuMultiStageSim` 的历史语义一致。历史为每个实际执行的步保存一行；活跃列表压缩在 E3 中明确推迟且尚未实现。

**后端模型**：CUDA 可用时，`Direct` 使用直接互感流水线。`Graph` 仅是图辅助：可捕获/重放互感 segment，而矩阵组装、求解器、力/状态编排和热处理仍在图外。`Persistent` 是一个请求；当前同步引擎没有常驻 control stream，因此它可能解析为 `Fallback`，多行 batch 也不能解析为 persistent。显式 `Fallback` 只使用 CPU，且不创建 CUDA context。

**示例——电容电压参数扫描**：

```cpp
#include <coilgun/coilgun_cuda.hpp>
#include <iostream>
#include <vector>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::SimBatch;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;
    using coilgun::simulation::TriggerConfig;
    using coilgun::simulation::TriggerMode;
    using coilgun::simulation::TerminationPolicy;

    // 共享几何
    DrivingCoil c1(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    int N = 5;
    SimBatch<EulerStepper> batch({c1, c2}, arm, N, 1e-6);

    // 配置每仿真激励
    double voltages[] = {300.0, 350.0, 400.0, 450.0, 500.0};
    for (int i = 0; i < N; ++i) {
        std::vector<std::unique_ptr<Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(voltages[i], 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(voltages[i], 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        batch.set_excitations(i, std::move(excs), triggers);
    }

    auto pol = TerminationPolicy::defaults();
    pol.max_steps = 500;  // 截断以加速扫描
    batch.run(pol);

    for (int i = 0; i < N; ++i) {
        auto& r = batch.result(i);
        std::cout << "V=" << voltages[i] << "V  →  v="
                  << r.summary.muzzle_velocity << " m/s\n";
    }
    return 0;
}
```

---

### GpuAdaptor（遗留/内部兼容）

```cpp
#include <coilgun/simulation/cuda/gpu_adaptor.hpp>

namespace coilgun::simulation::cuda {

struct CoilGeo { double ri, re, len, pos; int turns; };
struct FilGeo  { double ri, re, len; };

class GpuAdaptor {
public:
    GpuAdaptor();
    ~GpuAdaptor();
    GpuAdaptor(const GpuAdaptor&) = delete;
    GpuAdaptor& operator=(const GpuAdaptor&) = delete;
    GpuAdaptor(GpuAdaptor&&) noexcept;
    GpuAdaptor& operator=(GpuAdaptor&&) noexcept;

    void setup(const std::vector<DrivingCoil>& coils, const Armature& arm, int n_nodes);
    void setup_batch(const std::vector<DrivingCoil>& coils, const Armature& arm, int num_sims, int n_nodes);

    void upload_separation(const std::vector<double>& seps);
    void download_results(std::vector<double>& M_out, std::vector<double>& dM_out, int n_pairs);
    void upload_batch_separations(const std::vector<double>& seps);
    void download_batch_results(std::vector<double>& M_out, std::vector<double>& dM_out);

    const CoilGeo*  d_coils()   const;
    const FilGeo*   d_fils()    const;
    const double*   d_nodes()    const;
    const double*   d_weights()  const;
    double*         d_results_M();
    double*         d_results_dM();
    double*         d_seps();
    double*         d_batch_seps();
    double*         d_batch_results_M();
    double*         d_batch_results_dM();

    int n_stages() const; int n_fil() const; int n_nodes() const; int batch_size() const;
};

}
```

`GpuAdaptor` 是保留用于内部兼容的遗留、仅移动设备内存辅助类。当前的 `GpuSingleStageSim` 和 `GpuMultiStageSim` 都通过 `GpuEngine` 执行；只有 `SimBatch` 和引擎内部的遗留持久化后端实现仍使用 `GpuAdaptor`。大多数用户不应直接使用它。

下面的 `setup()`/上传契约只记录历史 adaptor 架构。它不是当前 `GpuEngine` 契约，也不保证当前单级引擎只上传一次不变数据或暴露这些设备缓冲区。

**`CoilGeo` / `FilGeo`** 是用于设备传输的紧凑 POD 结构体。它们分别镜像 `DrivingCoil` 和 `Armature` 的几何字段，展平后用于 GPU kernel 参数空间。

**遗留 `setup()`** 在历史 adaptor 架构中分配并将不变几何（线圈、丝元、GL 节点/权重）上传到设备内存，并在 per-step 操作前调用一次；这不是当前 `GpuEngine` 的保证。

**`setup_batch()`** 扩展 `setup()`，增加 `SimBatch` 所需的 per-sim 批量缓冲区。分配 `num_sims` 宽的分离和结果数组。

**`upload_separation()`** / **`download_results()`** 处理单个仿真每个步骤的主机↔设备数据移动。

**`upload_batch_separations()`** / **`download_batch_results()`** 处理批量仿真的数据移动。

设备指针访问器（`d_*()`）返回已分配的设备内存指针——由 CUDA kernel 直接使用。注意 `d_batch_*()` 指针仅在调用 `setup_batch()` 后有效。

**注意**：`GpuAdaptor` 是仅移动类型（删除拷贝）。析构函数通过 `cudaFree()` 释放所有设备分配。本节仅为遗留/内部兼容保留，不是当前单级 API 的权威说明。

---

### 内部设备头文件

`.cuh` 头文件定义 GPU kernel 使用的 `__host__ __device__` 函数。它们不面向用户代码。

```cpp
#include <coilgun/physics/elliptic.cuh>          // __host__ __device__ 椭圆积分
#include <coilgun/physics/mutual_inductance.cuh> // __host__ __device__ 丝级 M、dM/dz
```

这些头文仅能用 `nvcc` 编译（受 `#ifdef __CUDACC__` 保护）。它们提供了椭圆积分和丝级互感函数的 GPU 兼容内联实现。`coilgun_cuda.hpp` umbrella 头**不**包含它们——仅供内部使用。

---

### 遗留架构参考

下面的图和传输表描述历史 `GpuAdaptor` 路径，仅为兼容文档保留，不是当前 `GpuEngine` 的执行流程。当前单级和多级权威流程见上文：`Graph` 只捕获/重放 mutual-inductance segment，矩阵组装、Eigen 求解、力/状态更新和热处理仍在图外同步直接执行；`Fallback` 在 CPU 上执行完整物理步。

**每时间步的计算流程**：

```
┌─────────────────────────────────────────────┐
│ Host (CPU)                                   │
│   check_triggers() → extinguish_quiet()       │
│   填充映射分离值 / 门铃                         │
│   等待持久化 kernel（或启动逐对 kernel）       │
│   读取映射结果（或 cudaMemcpy D→H）            │
│   build_system_matrix [L - M_I]              │
│   Eigen LDLT 求解 → 新电流                    │
│   compute_force(F = Σ I_d × I_f × dM)        │
│   更新速度 / 位置                              │
│   更新电容电压 / 温度                          │
├─────────────────────────────────────────────┤
│ Device (GPU)                                  │
│   persistent_batch_kernel 或                 │
│   mutual_inductance_coil_pair_kernel          │
│     每 block: 512 threads × ~13 次循环        │
│     每次循环: 1 对椭圆积分                     │
│     shared memory 树形归约                     │
│     → 每对输出 1 个 double M, 1 个 double dM   │
└─────────────────────────────────────────────┘
```

**历史上构造时一次性上传的数据**（通过遗留 `GpuAdaptor`）：

| 缓冲区 | 大小 | 内容 |
|--------|------|---------|
| `d_coils_` | `n_stages × sizeof(CoilGeo)` | 每级 ri, re, length, position, turns |
| `d_fils_` | `N_fil × sizeof(FilGeo)` | 每丝元环 ri, re, length |
| `d_nodes_` | `9 × sizeof(double)` | Gauss-Legendre 积分节点 |
| `d_weights_` | `9 × sizeof(double)` | GL 积分权重 |

**历史上每步传输的数据**：

| 方向 | 大小 | 内容 |
|-----------|------|---------|
| H→D | `n_active × N_fil × sizeof(double)` | 电枢位置 → 分离值 |
| D→H | `n_stages × N_fil × 2 × sizeof(double)` | M1 和 dM1 矩阵 |

### 性能特征

| 规模 (S×F) | CPU (16核) | GPU (RTX 5080) | GPU 优势 |
|---|---|---|---|
| 1×10 | 19 s | 16 s | 1.2× |
| 2×10 | 58 s | 52 s | 1.1× |
| 25×45（典型） | ~5 min | ~2 min（估计） | ~2.5× |
| 50×200（高分辨率） | ~2 h | ~15 min（估计） | ~8× |

GPU 优势随问题规模增大而提升，因为 4D 积分 kernel（每对 6561 次椭圆积分求值）暴露大量并行性。小规模时 kernel launch overhead 和 PCIe 传输占主导。大规模时 GPU 计算吞吐量饱和。

### 线程安全

GPU 类是**单线程**的——不支持在同一实例上并发调用 `step()`。`SimBatch` 将主机侧仿真更新串行化；公开 GPU 仿真构造函数不接收 `cudaStream_t`，因此独立 stream 执行不属于此 API。

底层 CUDA kernel 对 host 是线程安全的：映射内存和 kernel launch 使用默认流，操作被串行化。

### 已知限制

| 限制 | 详情 |
|---|---|
| 自适应 GL 阶数（n_nodes=4/9） | 已移除。GPU 上使用 4 个 GL 节点会导致非确定性浮点漂移（B1）。 |
| 热模式 | `ThermalMode::Cpu` 使用 CPU 材料表更新；支持时 `ThermalMode::Gpu` 使用 GPU thermal workspace。解析后的模式和 `thermal_time_ms` 会记录在 `ExecutionReport` 中。 |
| 持久化 kernel | `GpuSingleStageSim`、`GpuMultiStageSim` 和 `SimBatch` 通过 `GpuEngine` 解析后端；只有在引擎能力和确定性策略允许时才尝试持久化执行，否则报告会标识 CPU 回退。 |
| CUDA Graphs | 仅实现 mutual-inductance segment 的捕获。矩阵组装、Eigen 求解、力/状态更新和热处理仍在图外执行；报告会标识实际 Graph 路径，但不宣称完整流水线已被捕获。 |
| 双缓冲 | 未实现。持久化结果在 CPU LDLT 求解前同步，不承诺 CPU/GPU 重叠执行。 |

测量目标 `bench_gpu_engine` 同时包含 CPU Reference 基线，并记录依赖机器的
墙钟、求解器、热路径、传输和 Graph 捕获观测值。结果见
`docs/benchmarks/2026-07-19-unified-gpu-engine.md`。

---

## 构建系统参考

### CMake 选项

| 选项 | 默认值 | 说明 |
|--------|---------|-------------|
| `COILGUN_BUILD_TESTS` | `ON` | 构建单元测试和集成测试 |
| `COILGUN_BUILD_GENERATOR` | `OFF` | 构建 T(q,p) 查表生成工具 |
| `COILGUN_ENABLE_CUDA` | `OFF` | 构建 GPU 加速后端 (`libcoilgun_cuda.a`) |

### CMake 预设

| 预设 | 生成器 | 构建类型 | 标志 | 备注 |
|--------|-----------|------------|-------|-------|
| `ninja-debug` | Ninja | Debug | `-march=native` | 启用测试，生成 compile_commands.json |
| `ninja-release` | Ninja | Release | `-march=native -O3` | — |
| `make-debug` | Unix Makefiles | Debug | `-march=native` | 启用测试，生成 compile_commands.json |
| `ninja-cuda-debug` | Ninja | Debug | `-march=native` | CUDA 启用，测试启用 |

### 测试预设

```sh
ctest --preset debug   # 运行全部 18 个已配置测试
```
