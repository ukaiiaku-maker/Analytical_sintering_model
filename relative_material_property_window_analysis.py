#!/usr/bin/env python3
"""Exact-first attribution and threshold extraction for material windows."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

OUT=Path("results/relative_material_property_window_attribution");SRC=OUT/"source_tables";RAW=OUT/"raw_outputs"
FEATURES=["Q_nuc_delta_kJ","Q_exchange_delta_kJ","Q_transport_delta_kJ","Q_growth_delta_kJ","Q_PR_delta_kJ","Q_closed_delta_kJ",
          "k_nuc_factor","k_exchange_factor","k_transport_factor","k_growth_factor","k_PR_factor","k_closed_factor",
          "Theta_nuc","f_nuc","f_exchange","f_transport","I_low_slow","I_low_PR_slow","Pi_PR","S_closed_growth","M_PR_closed",
          "Q_nuc_minus_Q_growth_kJ","Q_nuc_minus_Q_PR_kJ","Q_nuc_minus_Q_transport_kJ","Q_closed_minus_Q_growth_kJ","Q_PR_minus_Q_closed_kJ",
          "log10_kclosed_over_kgrowth","log10_kPR_over_kgrowth","log10_knuc_over_ktransport"]


def bools(x):return x.astype("boolean").fillna(False).to_numpy(bool)


def main():
    screen=pd.read_csv(SRC/"material_property_window_scorecard.csv",low_memory=False)
    exact=pd.read_csv(SRC/"material_property_window_exact_promotions.csv")
    exact=exact.merge(screen,on="property_id",how="left",suffixes=("","_screenrow"))
    fp=bools(exact.fast_firing_pass_exact);tp=bools(exact.two_step_pass_exact)
    exact["both_pass_exact"]=fp&tp
    # Sensitivity rankings use exact outcomes when observed, with screen rows
    # retained as explicitly labeled feasibility evidence.
    targets=(("fast_firing_pass_exact",fp.astype(float),"exact"),("two_step_pass_exact",tp.astype(float),"exact"),
             ("both_pass_exact",(fp&tp).astype(float),"exact"),("R_fast_exact",exact.R_fast_exact,"exact"),
             ("reduction_TS_exact",exact.reduction_TS_exact,"exact"),("Chen_window_width_C_exact",exact.Chen_window_width_C_exact,"exact"),
             ("fast_firing_pass_screen",screen.fast_firing_pass_screen.astype(float),"screen"),
             ("two_step_pass_screen",screen.two_step_pass_screen.astype(float),"screen"))
    ranks=[]
    for target,y,level in targets:
        for feature in FEATURES:
            frame=exact if level=="exact" else screen
            yy=pd.Series(y,index=frame.index);mask=frame[feature].notna()&yy.notna()
            corr=float(frame.loc[mask,feature].rank().corr(yy[mask].rank())) if mask.sum()>3 and frame.loc[mask,feature].nunique()>1 else np.nan
            ranks.append(dict(target=target,feature=feature,evidence_level=level,n=int(mask.sum()),spearman_r=corr,absolute_rank_score=abs(corr) if np.isfinite(corr) else np.nan))
    rankings=pd.DataFrame(ranks).sort_values(["target","absolute_rank_score"],ascending=[True,False]);rankings.to_csv(SRC/"material_property_sensitivity_rankings.csv",index=False)

    groups={"fast_pass_exact":fp,"two_step_pass_exact":tp,"both_pass_exact":fp&tp}
    threshold_features=["Q_nuc_minus_Q_growth_kJ","Q_nuc_minus_Q_PR_kJ","Q_nuc_minus_Q_transport_kJ","Q_closed_minus_Q_growth_kJ",
                        "Q_PR_minus_Q_closed_kJ","log10_kclosed_over_kgrowth","log10_kPR_over_kgrowth","log10_knuc_over_ktransport",
                        "Theta_nuc","f_nuc","I_low_PR_slow","S_closed_growth","M_PR_closed","A_closed_fraction"]
    thresholds=[]
    for group,mask in groups.items():
        q=exact[mask]
        for feature in threshold_features:
            thresholds.append(dict(group=group,feature=feature,n=len(q),minimum=q[feature].min(),maximum=q[feature].max(),median=q[feature].median(),
                                   evidence="exact promoted rows; coverage-limited, not universal bounds"))
    pd.DataFrame(thresholds).to_csv(SRC/"dimensionless_thresholds.csv",index=False)

    activation=[]
    for group,mask in groups.items():
        q=exact[mask]
        for feature in ("Q_nuc_delta_kJ","Q_exchange_delta_kJ","Q_transport_delta_kJ","Q_growth_delta_kJ","Q_PR_delta_kJ","Q_closed_delta_kJ"):
            activation.append(dict(group=group,parameter=feature,n=len(q),minimum_delta_kJ_mol=q[feature].min(),maximum_delta_kJ_mol=q[feature].max(),
                                   total_observed_width_kJ_mol=q[feature].max()-q[feature].min(),evidence="exact promotion envelope"))
    activation=pd.DataFrame(activation);activation.to_csv(SRC/"relative_activation_energy_window.csv",index=False)

    ingredients=[
      ("nucleation-limited waiting","fast firing","necessary","nucleation-facile ablation removes exact fast cases"),
      ("PR redistribution","fast firing","non-causal in current envelope","PR-off may preserve fast firing"),
      ("exchange and transport completion","fast firing","necessary in combination","both paths must attain matched-density interval"),
      ("PR-prepared closed store","two-step","necessary for candidate 693168","no-PR damage destroys frozen candidate"),
      ("closed shrinkage","two-step","necessary","supports density above 0.95 and lower boundary"),
      ("finite closed accommodation","two-step","necessary and calibration-sensitive","infinite-accommodation ablation destroys joint result"),
      ("thermally activated migration","two-step","necessary for upper boundary","growth failure brackets success band"),
      ("high-density attainment","both","shared constraint","unattained intervals never score"),
      ("large grain-size separation","two-step","outcome, not artifact","allowed when attained and numerically supported"),
    ]
    pd.DataFrame(ingredients,columns=["ingredient","behavior","attribution","evidence"]).to_csv(SRC/"mechanism_attribution_summary.csv",index=False)

    # Reduced family transfer: apply the exact-693168 OAT survival fractions to
    # each frozen Tier-B margin.  It is explicitly not an exact family sweep.
    fam=pd.read_csv("results/reframe_tierB_experimental_plausibility/tierB_candidate_reinterpretation.csv")
    oat=exact[exact.design_stage=="OAT"]
    base_red=float(fam.loc[fam.candidate_id==693168,"median_reduction"].iloc[0]);base_w=float(fam.loc[fam.candidate_id==693168,"window_width_C"].iloc[0])
    family=[]
    for _,r in fam.iterrows():
        margin=min(float(r.median_reduction)/base_red,float(r.window_width_C)/base_w)
        family.append(dict(candidate_id=int(r.candidate_id),base_reduction=r.median_reduction,base_window_width_C=r.window_width_C,
                           closed_fraction_at_switch=r.closed_fraction_at_switch,robustness_count=r.robustness_count,
                           transferred_material_margin=margin,predicted_Qclosed_window_kJ_mol=200*min(margin,1),
                           predicted_Qgrowth_window_kJ_mol=200*min(margin,1),predicted_prefactor_decades=3*min(margin,1),
                           reduced_family_class="robust comparator" if margin>=.5 and r.robustness_count>0 else "calibration-sensitive",
                           evidence_level="reduced transfer from exact 693168 OAT; topology frozen; not exact family validation"))
    pd.DataFrame(family).to_csv(SRC/"tierB_family_material_window_comparison.csv",index=False)

    counts=pd.DataFrame([dict(screened_rows=len(screen),exact_union_rows=len(exact),exact_fast_rows=int(exact.fast_firing_pass_exact.notna().sum()),
                                  exact_two_step_rows=int(exact.two_step_pass_exact.notna().sum()),fast_only=int((fp&~tp).sum()),two_step_only=int((tp&~fp).sum()),
                                  both_pass=int((fp&tp).sum()),neither=int((~fp&~tp).sum()),surrogate_both=int((screen.classification_screen=="both_pass").sum()))])
    counts.to_csv(SRC/"material_property_window_classification_summary.csv",index=False)
    state=dict(status="complete",**counts.iloc[0].to_dict(),screen_is_exact=False,model_layers_coupled=False)
    (OUT/"analysis_run_state.json").write_text(json.dumps(state,indent=2)+"\n")
    print(json.dumps(state,indent=2))


if __name__=="__main__":main()
