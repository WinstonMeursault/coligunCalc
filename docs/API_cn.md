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

    extern const MaterialProperties COPPER;        // 1.75e-8 Ω·m, 0.0039 K⁻¹, 8960 kg/m³
    extern const MaterialProperties ALUMINUM;      // 2.82e-8 Ω·m, 0.0040 K⁻¹, 2700 kg/m³

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
L = μ₀ · nc² · ri³ · T(q, p)
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
// 查表（q=5, p=3 均在范围内）：
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
    double force;                         // 电磁力
    std::vector<double> filament_temperatures; // 各丝温度（仅启用热耦合时）
};

struct SimSummary {
    double muzzle_velocity;    // 出口速度 (m/s)
    double total_time;         // 总耗时 (s)
    double max_force;          // 最大力 (N)
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

template<typename SP = EulerStepper>
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
    Reference   = 0,   // 全部参考性计算：自感精确计算、无距离截断、固定 9 点 GL 积分
    LookupTable = 1,   // 仅启用 T(q,p) 查表加速自感计算，其他为参考性计算
    Full        = 2    // 全优化：查表 + 距离截断 + 自适应 GL 阶数（近 9 远 4）
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

### 多级仿真结果

```cpp
#include <coilgun/simulation/multi_stage_result.hpp>

struct StepSnapshot {
    double time;                    // 仿真时间 s
    double arm_position;            // 电枢中心位置 m
    double arm_velocity;            // 电枢速度 m/s
    double force;                   // 轴向力 N
    std::vector<double> filament_currents;       // 各纤维电流 A
    std::vector<double> filament_temperatures;   // 各纤维温度 K（仅热模式）
};

struct MultiStageStep {
    StepSnapshot          state;
    std::vector<double>   cap_voltages;     // [n_stages]，未触发级为 0
    std::vector<double>   coil_currents;    // [n_stages]，未触发级为 0
};

struct PerStageSummary {
    int    stage_index;           // 零基级号
    double trigger_time;          // 触发时刻 s
    double trigger_position;      // 触发时电枢位置 m
    double peak_current;          // 该级峰值电流 A
    double max_force;             // 该级贡献的最大力 N
    double energy_depleted;       // 电容消耗能量 J
    int    step_count_active;     // 该级活跃步数
};

struct MultiStageSummary {
    double                        muzzle_velocity;      // 出口速度 m/s
    double                        total_time;           // 仿真总时间 s
    double                        max_force;            // 全局峰值力 N
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

template<typename SP = EulerStepper>
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

**优化**：`opt_level` 参数控制运行时计算优化，独立于组件构造时的自感计算模式（后者在构造 DrivingCoil/Armature 时已确定）。在 `Reference` 级别，所有互感计算使用固定 9 点 Gauss-Legendre 求积且无距离截断。在 `Full` 级别，远处线圈-纤维对被跳过，近处使用 9 点、远处使用 4 点求积。

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
    bool    enable_profiling  = false; ///< 插入 NVTX 范围标记，用于 nsight 分析。
};

}
```

| 字段 | 说明 | 默认值 |
|-------|---------|---------|
| `device_id` | 选择目标物理 GPU（多 GPU 系统相关）。 | `0` |
| `threads_per_block` | 4D GL 积分 kernel 的每 block 线程数。必须是 2 的幂且 ≤ 1024。较大值减少每线程循环次数但增加共享内存压力。 | `512` |
| `max_batch_sims` | `SimBatch` 中的最大仿真数。缓冲区预分配到此大小。超出则抛出 `std::runtime_error`。 | `256` |
| `enable_profiling` | 为 true 时，在每个 `step()` 周围插入 NVTX push/pop 范围用于 GPU 时间线分析。 | `false` |

**示例**：

```cpp
coilgun::simulation::cuda::GpuBackend be;
be.device_id = 0;
be.threads_per_block = 1024;  // 最大化并行度
```

---

### GpuOptLevel

```cpp
#include <coilgun/simulation/cuda/gpu_backend.hpp>

namespace coilgun::simulation::cuda {

enum class GpuOptLevel {
    Standard = 0,   ///< n_nodes=9，无距离截断。用于调试/验证。
    Full     = 1,   ///< 距离截断（>10× 线圈长度）+ 固定 n_nodes=9。生产环境。
};

}
```

| 等级 | 距离截断 | GL 阶数 | 行为 |
|:-----:|:---------------:|:--------:|-----------|
| `Standard` | 否 | n_nodes = 9 | 所有 (stage, filament) 对都处理，即使 stage 距离电枢很远。kernel 使用 9 个 GL 节点（6561 个积分点）。最安全的设置。 |
| `Full` | 是 | n_nodes = 9 | 距离电枢超过 10× 线圈长度的 stage 被完全跳过。在典型距离下所有活跃 stage 都在范围内。 |

**注意**：GPU 后端不支持自适应 GL 阶数（n_nodes=4/9 混合）。GPU 上使用 n_nodes=4 会导致共享内存归约的非确定性浮点漂移（仅 8 个 warp 中的 3 个有实际工作）。CPU `OptimizationLevel::Full` 同时使用距离截断和自适应积分；GPU `GpuOptLevel::Full` 仅使用距离截断。

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
        const GpuBackend&                backend = {}
    );

    const SimStep&   step();
    const SimResult& run();
    const SimResult& run(const TerminationPolicy& p);
    void             reset();

    const SimResult& result()     const;
    const SimState&  state()      const;
    double  dt()         const;
    int     step_count() const;
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

**方法**：

| 方法 | 行为 |
|--------|-----------|
| `step()` | 推进一个时间步。对每个丝元对发射 CUDA kernel，通过 Eigen LDLT 求解 ODE，更新运动学，记录诊断信息。返回新记录的 `SimStep`。 |
| `run()` | 运行至默认终止条件。阻塞直到激励结束（续流二极管导通且电流衰减完毕）。 |
| `run(policy)` | 运行至自定义 `TerminationPolicy`。 |
| `reset()` | 恢复至初始状态。所有电流归零，电枢位置/速度恢复至构造时值，重置激励，清空结果历史。 |
| `result()` | 最后一次 `run()` 之后的结果。包含 `history`（SimStep 向量）和 `summary`（SimSummary）。 |
| `state()` | 当前内部状态——实时引用。 |
| `dt()` | 固定时间步长。 |
| `step_count()` | 自构造或上次 `reset()` 以来的步数。 |

**内部实现**：类内部持有一个 `GpuAdaptor`，在构造时将线圈几何、丝元离散化和 GL 积分节点一次性上传到设备内存。每次 `step()` 调用将电枢位置传至设备，为每个（线圈，丝元）对发射 4D 积分 kernel，并回读 `M1`/`dM1` 矩阵。

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

公开 API 与 `SingleStageSim` 完全一致。从 CPU 迁移到 GPU 仅需修改类名和 include 路径。

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
        const GpuBackend&                         backend = {}
    );

    const MultiStageStep&  step();
    const MultiStageResult& run();
    const MultiStageResult& run(const TerminationPolicy& p);
    void                    reset();

    const MultiStageResult& result()      const;
    const MultiStageState&  state()       const;
    double  dt()              const;
    int     step_count()      const;
    int     num_stages()      const;
};

} // namespace
```

**构造约束**：`coils.size() == excitations.size()`，`trigger_configs.size() == coils.size() - 1`，`coils.size() ≤ kMaxStages (50)`。违反任一则抛出 `std::invalid_argument`。

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
| `backend` | GPU 后端配置 | `{}`（默认值） |

**方法**——签名与 CPU 版 `MultiStageSim` 完全一致：

| 方法 | 行为 |
|--------|-----------|
| `step()` | 推进一个时间步。检查触发器、静默已结束的 stage、为所有活跃的 (stage, filament) 对发射 CUDA kernel、求解 ODE、更新运动学。 |
| `run()` / `run(policy)` | 运行至终止条件。 |
| `reset()` | 恢复所有状态至初始条件。 |
| `result()` / `state()` / `dt()` / `step_count()` / `num_stages()` | 查询。 |

**内部实现**：使用 CUDA kernel 计算所有活跃 stage-丝元对的 M1/dM1 矩阵。丝元间互感矩阵 `[M]` 在构造时由 CPU 预计算（仿真过程中不变）。线圈间互感 `M_cc` 同样在构造时由 CPU 预计算，填充系统矩阵的 coil-coil 非对角块。LDLT 线性求解、力计算和运动学更新保持在 CPU 上。ODE 维度为 `n_stages + N_filaments`。未激活的 stage 的行/列设为单位矩阵。

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
| 优化等级 | `OptimizationLevel`（3 级） | `GpuOptLevel`（2 级） |
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
        const GpuBackend&                     backend = {});

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
};

} // namespace
```

`SimBatch` 是**参数扫描**的容器——运行共享相同线圈和电枢几何但激励参数（电压、电容）和/或触发位置不同的多组仿真。

**约束**：所有仿真必须共享**完全相同**的线圈几何和丝元离散化（`m × n`）。传入构造函数的 `coils` 和 `armature` 为所有仿真共享。每个仿真的激励源和触发配置通过 `set_excitations()` 设置。

**构造参数**：

| 参数 | 说明 |
|-----------|-------------|
| `coils` | 共享驱动线圈几何（拷贝）。所有仿真相同。 |
| `armature` | 共享电枢几何和丝元离散化（拷贝）。 |
| `num_sims` | 并发仿真数。必须 ≤ GpuBackend 中的 `max_batch_sims`。 |
| `dt` | 固定时间步长（s），所有仿真共享。 |
| `backend` | GPU 后端配置。 |

**方法**：

| 方法 | 行为 |
|--------|-----------|
| `set_excitations(sim_id, excitations, triggers)` | 为仿真 `sim_id` 配置激励源和触发设置。必须在 `run()` 之前对每个仿真调用。 |
| `run()` / `run(policy)` | 同时运行所有仿真。每个仿真独立运行——最先到达终止条件的仿真不阻塞其他仿真。 |
| `result(sim_id)` | 获取特定仿真的 `MultiStageResult`。 |
| `num_sims()` | 本批次中的仿真数。 |

**执行模型**：`SimBatch` 内部管理一个 `GpuAdaptor` 和一个 per-sim 状态对象数组。每个时间步先对所有仿真检查触发器和静默 stage，然后通过 GPU kernel launch 循环计算所有活跃对的 M1/dM1，最后在 CPU 上对每个仿真独立求解 ODE 和更新运动学。

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

**注意**：当前实现对每个 step 的每个 (stage, filament) 对、每个仿真发射一个 CUDA kernel（对典型规模 launch overhead ~100 µs/步）。功能上等价于用循环调用 `GpuMultiStageSim`，但提供了统一的 API 和共享几何。持久化 kernel 优化（每步零 kernel launch）已设计，计划在后续版本中实现。

---

### GpuAdaptor（高级）

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

`GpuAdaptor` 管理 GPU 加速仿真的设备内存。被 `GpuSingleStageSim`、`GpuMultiStageSim` 和 `SimBatch` 内部使用。大多数用户不需要直接与之交互。

**`CoilGeo` / `FilGeo`** 是用于设备传输的紧凑 POD 结构体。它们分别镜像 `DrivingCoil` 和 `Armature` 的几何字段，展平后用于 GPU kernel 参数空间。

**`setup()`** 分配并将不变几何（线圈、丝元、GL 节点/权重）上传到设备内存。必须在任何 per-step 操作之前调用一次。

**`setup_batch()`** 扩展 `setup()`，增加 `SimBatch` 所需的 per-sim 批量缓冲区。分配 `num_sims` 宽的分离和结果数组。

**`upload_separation()`** / **`download_results()`** 处理单个仿真每个步骤的主机↔设备数据移动。

**`upload_batch_separations()`** / **`download_batch_results()`** 处理批量仿真的数据移动。

设备指针访问器（`d_*()`）返回已分配的设备内存指针——由 CUDA kernel 直接使用。注意 `d_batch_*()` 指针仅在调用 `setup_batch()` 后有效。

**注意**：`GpuAdaptor` 是仅移动类型（删除拷贝）。析构函数通过 `cudaFree()` 释放所有设备分配。

---

### 内部设备头文件

`.cuh` 头文件定义 GPU kernel 使用的 `__host__ __device__` 函数。它们不面向用户代码。

```cpp
#include <coilgun/physics/elliptic.cuh>          // __host__ __device__ 椭圆积分
#include <coilgun/physics/mutual_inductance.cuh> // __host__ __device__ 丝级 M、dM/dz
```

这些头文仅能用 `nvcc` 编译（受 `#ifdef __CUDACC__` 保护）。它们提供了椭圆积分和丝级互感函数的 GPU 兼容内联实现。`coilgun_cuda.hpp` umbrella 头**不**包含它们——仅供内部使用。

---

### 架构

**每时间步的计算流程**：

```
┌─────────────────────────────────────────────┐
│ Host (CPU)                                   │
│   check_triggers() → extinguish_quiet()       │
│   对每个 (stage, filament) 对:               │
│     计算分离 → CUDA kernel launch            │
│   cudaDeviceSynchronize()                    │
│   cudaMemcpy D→H: M1_mat, dM1_mat            │
│   build_system_matrix [L - M_I]              │
│   Eigen LDLT 求解 → 新电流                    │
│   compute_force(F = Σ I_d × I_f × dM)        │
│   更新速度 / 位置                              │
│   更新电容电压 / 温度                          │
├─────────────────────────────────────────────┤
│ Device (GPU)                                  │
│   mutual_inductance_coil_pair_kernel          │
│     每 block: 512 threads × ~13 次循环        │
│     每次循环: 1 对椭圆积分                     │
│     shared memory 树形归约                     │
│     → 每对输出 1 个 double M, 1 个 double dM   │
└─────────────────────────────────────────────┘
```

**构造时一次性上传的数据**（通过 `GpuAdaptor`）：

| 缓冲区 | 大小 | 内容 |
|--------|------|---------|
| `d_coils_` | `n_stages × sizeof(CoilGeo)` | 每级 ri, re, length, position, turns |
| `d_fils_` | `N_fil × sizeof(FilGeo)` | 每丝元环 ri, re, length |
| `d_nodes_` | `9 × sizeof(double)` | Gauss-Legendre 积分节点 |
| `d_weights_` | `9 × sizeof(double)` | GL 积分权重 |

**每步传输的数据**：

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

GPU 类是**单线程**的——不支持在同一实例上并发调用 `step()`。多仿真并发通过 `SimBatch`（串行化）或创建独立的 `GpuMultiStageSim` 实例并绑定到不同的 `cudaStream_t` 来实现。

底层 CUDA kernel 对 host 是线程安全的：映射内存和 kernel launch 使用默认流，操作被串行化。

### 已知限制

| 限制 | 详情 |
|---|---|
| 自适应 GL 阶数（n_nodes=4/9） | 已移除。GPU 上使用 4 个 GL 节点会导致非确定性浮点漂移（B1）。 |
| 热模式 | 温度更新在 CPU 上运行（单线程）。该计算 <1% 运行时间。 |
| Batch kernel | 真正的 batch kernel（单次 grid launch 处理所有仿真）已设计（P0）但被 nvcc 跨模块 struct-layout 问题阻塞。持久化 kernel 优化（D1-D2）是计划中的替代方案。 |
| CUDA Graphs | 未实现。需要每步固定 n_active（无动态触发/终止分支）。 |

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
ctest --preset debug   # 运行全部 17 套测试
```
