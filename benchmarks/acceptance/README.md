# Acceptance 独立验收数据集（2026-08-18 建立）

6 个独立靶点，与 AB-Bind core_v1 / 现有样本外集（1DVF/1AK4/1KTZ/3K2M）零重叠。
结构：RCSB 原件（structures/）；突变与实验 ddG：SKEMPI 2.0（RT·ln(Kd_mut/Kd_wt)，
行级 assay 温度；同突变多测量取均值）。

| 靶点 | 复合物 | 分辨率 | V1（电荷守恒） | V2.1a（电荷变化） |
|---|---|---|---|---|
| 4I77 | Lebrikizumab Fab–IL-13 | 1.9 | 11 | 1 |
| 2BDN | 11K2 scFv–MCP-1 | 2.53 | 2 | 10 |
| 1NMB | NC10 Fab–流感 N9 NA | 2.2 | 6 | 1 |
| 1AHW | 5G9 Fab–组织因子 | 3.0 | 8 | 1 |
| 1IAR | IL-4–IL-4Rα | 2.3 | 11 | 25 |
| 1OGA | JM22 TCR–HLA-A2/Tax | 1.4 | 31 | 17 |

- 结构 QC：声明链全部在场、突变位点 WT 身份零错配（2026-08-18 验证）。
- 逐靶点 mutations_v1.csv / mutations_v21a.csv + 对应 experimental_ddg 文件。
- 抗原多样性：细胞因子×2、趋化因子、病毒抗原、凝血因子受体、pMHC/TCR。
