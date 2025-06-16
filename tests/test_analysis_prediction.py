# -*- coding: utf-8 -*-
"""
Created on Mon May  6 16:11:38 2024

@author: veenstra
"""

import pytest
import os
import numpy as np
import pandas as pd
import hatyan
import pytz
import datetime as dt
from hatyan.analysis_prediction import MatrixConditionTooHigh

dir_tests = os.path.dirname(__file__) #F9 doesnt work, only F5 (F5 also only method to reload external definition scripts)
dir_testdata = os.path.join(dir_tests,'data_unitsystemtests')


@pytest.mark.unittest
def test_analysis_settings():
    
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1]]
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    
    ts_comp_nfac1_fualltimes0_xfac1_peryear0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, analysis_perperiod=False)
    
    ts_comp_onecomp = hatyan.analysis(ts=ts_measurements_group0, const_list=['M2'], nodalfactors=True, fu_alltimes=True, xfac=True)
    
    ts_comp_nfac1_fualltimes1_xfac1 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True)
    ts_comp_nfac1_fualltimes0_xfac1 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True)
    ts_comp_nfac1_fualltimes1_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=False)
    ts_comp_nfac1_fualltimes0_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=False)
    ts_comp_nfac0_fualltimes0_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=False, fu_alltimes=False, xfac=False)
    
    comp_A_expected = np.array([9.92687566e-04, 3.11039747e-02, 9.73652407e-02, 6.68645368e-02,
                                4.02714135e-02, 2.44247279e-02, 1.34401343e-01, 2.69848016e-01,
                                1.74658516e+00, 1.33037114e-01, 4.78882284e-01, 4.42634321e-02,
                                1.93477077e-02, 4.29425377e-02, 1.28103067e-01, 8.89165480e-02,
                                4.66780541e-02, 8.50788597e-02, 9.01287840e-02, 3.06546626e-02,
                                4.72786278e-02, 1.47588322e-02])
    
    comp_phi_expected = np.array([  0.        , 141.78046951, 188.22394202,   7.35166063,
                                  279.19176166, 122.65613222, 159.77053294,  33.54516691,
                                   59.02752661, 256.04742485, 117.45538945, 345.62388999,
                                  200.90891452,  92.80114988, 115.31391482, 176.58434752,
                                   77.32631928, 103.4830287 , 153.94793646, 110.76924584,
                                  157.3933541 , 220.41306556])
    
    assert (np.abs(ts_comp_nfac1_fualltimes0_xfac1_peryear0-ts_comp_nfac1_fualltimes0_xfac1) < 10E-9).all().all()
    assert (np.abs(ts_comp_nfac1_fualltimes0_xfac1_peryear0['phi_deg'].values-comp_phi_expected) < 10E-9).all()
    assert (np.abs(ts_comp_nfac1_fualltimes0_xfac1_peryear0['A'].values-comp_A_expected) < 10E-9).all()

    assert np.allclose(ts_comp_nfac1_fualltimes0_xfac1_peryear0.loc["M2"].values, np.array([ 1.74658516, 59.02752661]))
    assert np.allclose(ts_comp_onecomp.loc["M2"].values, np.array([ 1.74229247, 58.77828766]))
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac1.loc["M2"].values, np.array([ 1.74649803, 59.01827634]))
    assert np.allclose(ts_comp_nfac1_fualltimes0_xfac1.loc["M2"].values, np.array([ 1.74658516, 59.02752661]))
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac0.loc["M2"].values, np.array([ 1.7623892 , 59.01684785]))
    assert np.allclose(ts_comp_nfac1_fualltimes0_xfac0.loc["M2"].values, np.array([ 1.76257914, 59.02752661]))
    assert np.allclose(ts_comp_nfac0_fualltimes0_xfac0.loc["M2"].values, np.array([ 1.72889408, 57.192191  ]))


@pytest.mark.unittest
def test_analysis_settings_extended():
    
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1,2]]
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    
    ts_comp_nfac1_fualltimes1_xfac1_peryear0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, analysis_perperiod=False)
    
    ts_comp_nfac1_fualltimes1_xfac1_permonth0, comp_allperiods = hatyan.analysis(ts=ts_measurements_group0, const_list='month', 
                                                                                 nodalfactors=True, fu_alltimes=True, xfac=True, 
                                                                                 analysis_perperiod='M', return_allperiods=True)
    assert comp_allperiods.shape == (22, 48)
    
    ts_comp_nfac1_fualltimes1_xfac1 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, analysis_perperiod='Y')
    ts_comp_nfac1_fualltimes0_xfac1 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, analysis_perperiod='Y')
    ts_comp_nfac1_fualltimes1_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=False, analysis_perperiod='Y')
    ts_comp_nfac1_fualltimes0_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=False, analysis_perperiod='Y')
    ts_comp_nfac0_fualltimes0_xfac0 = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=False, fu_alltimes=False, xfac=False, analysis_perperiod='Y')
    
    # assert if the results are somewhat close to each other
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac1_peryear0, ts_comp_nfac1_fualltimes1_xfac1_permonth0, rtol=1e-1, atol=1e-1)
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac1, ts_comp_nfac1_fualltimes0_xfac1, rtol=1e-2, atol=1e-3)
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac1, ts_comp_nfac1_fualltimes0_xfac1, rtol=1e-2, atol=1e-2)
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac1, ts_comp_nfac1_fualltimes1_xfac0, rtol=1e-2, atol=1e-2)
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac1, ts_comp_nfac1_fualltimes0_xfac0, rtol=1e-2, atol=1e-2)
    assert np.allclose(ts_comp_nfac1_fualltimes1_xfac1, ts_comp_nfac0_fualltimes0_xfac0, rtol=1e1, atol=1e0)


@pytest.mark.unittest
def test_analysis_settings_perperiod():
    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0)
    
    comp_mean, comp_all = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, 
                                          analysis_perperiod="M", return_allperiods=True)
    assert len(comp_mean.attrs) > 0
    assert len(comp_all.attrs) > 0


@pytest.mark.unittest
def test_analysis_component_splitting():
    current_station = 'D15'
    file_ts = os.path.join(dir_testdata,'%s_obs1.txt'%(current_station))
    cs_comps = pd.DataFrame({'CS_comps_from':['K1','N2','2MN2','S2','S2'],
                             'CS_comps_derive':['P1','NU2','LABDA2','K2','T2'],
                             'CS_ampfacs':[0.33,0.22,0.48,0.29,0.05],
                             'CS_degincrs':[-11,-24,174,1,-24]})
    
    ts_meas = hatyan.read_dia(filename=file_ts, station=current_station)
    comp_cs_atonce = hatyan.analysis(ts=ts_meas, const_list='month', 
                                     analysis_perperiod=False, cs_comps=cs_comps)
    comp_cs_peryear = hatyan.analysis(ts=ts_meas, const_list='month', 
                                      analysis_perperiod="Y", cs_comps=cs_comps)
    comp_cs_permonth = hatyan.analysis(ts=ts_meas, const_list='month', 
                                       analysis_perperiod="M", cs_comps=cs_comps)
    
    assert np.allclose(comp_cs_atonce.loc["M2"].values, [  0.6222659 , 187.24099172])
    assert np.allclose(comp_cs_peryear.loc["M2"].values, [  0.6222659 , 187.24099172])
    assert np.allclose(comp_cs_permonth.loc["M2"].values, [  0.6222659 , 187.24099172])


@pytest.mark.unittest
def test_analysis_foreman():
    current_station = 'VLISSGN'
    
    #comp0
    file_data_comp0 = [os.path.join(dir_testdata,'%s_obs%i.txt'%(current_station, file_id)) for file_id in [1]]
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0, station=current_station)
    
    comp_frommeas_schu = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, source='schureman')
    comp_frommeas_for = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=True, xfac=True, source='foreman')
    schu_comp_phi_expected =  np.array([  0.        , 141.68386391, 188.27101759,   7.3992028 ,
                                       279.17807514, 122.81257343, 159.77061119,  33.51788732,
                                        59.01827634, 256.05796043, 117.47528038, 345.63438447,
                                       201.05307505,  92.81986396, 115.2820895 , 176.56952991,
                                        77.30925229, 103.43837741, 153.90952097, 110.64017054,
                                       157.30038107, 220.24021917])
    schu_comp_A_expected =    np.array([9.93677030e-04, 3.11102518e-02, 9.74940999e-02, 6.69089179e-02,
                                       4.01204828e-02, 2.43192049e-02, 1.34386805e-01, 2.69878623e-01,
                                       1.74649803e+00, 1.32976762e-01, 4.79019453e-01, 4.42086557e-02,
                                       1.92695651e-02, 4.29786476e-02, 1.28078709e-01, 8.89513755e-02,
                                       4.67055728e-02, 8.50832246e-02, 9.01742533e-02, 3.06526796e-02,
                                       4.72907179e-02, 1.47366661e-02])
    for_comp_phi_expected =  np.array([  0.        , 321.7403689 ,   8.15729584, 187.40999172,
                                       279.42669414, 122.65017113, 160.05967765,  33.32750808,
                                        59.01572119, 256.18495318, 117.35826025, 345.41402817,
                                       201.17575979,  92.60563462, 115.27544986, 176.43679913,
                                        77.08771066, 103.42351654, 153.76456299, 110.600143  ,
                                       157.16101443, 220.13261581])
    for_comp_A_expected =    np.array([9.96058525e-04, 3.01459094e-02, 9.72007804e-02, 6.69287134e-02,
                                       3.99816712e-02, 2.42471615e-02, 1.33778086e-01, 2.75247289e-01,
                                       1.76155491e+00, 1.39007149e-01, 4.85825897e-01, 4.40858825e-02,
                                       1.92198515e-02, 4.29817009e-02, 1.29452749e-01, 9.04975848e-02,
                                       4.66816274e-02, 8.62181973e-02, 9.27486375e-02, 3.13117938e-02,
                                       4.82570962e-02, 1.46874202e-02])
    
    assert (np.abs(comp_frommeas_schu['phi_deg'].values-schu_comp_phi_expected) < 10E-9).all()
    assert (np.abs(comp_frommeas_schu['A'].values-schu_comp_A_expected) < 10E-9).all()
    assert (np.abs(comp_frommeas_for['phi_deg'].values-for_comp_phi_expected) < 10E-9).all()
    assert (np.abs(comp_frommeas_for['A'].values-for_comp_A_expected) < 10E-9).all()


@pytest.mark.unittest
def test_analysis_settings_timeseries_tooshort():
    file_dia = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_pd = hatyan.read_dia(file_dia)
    
    with pytest.raises(MatrixConditionTooHigh) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list="all_schureman_originalorder")
    assert "condition of xTx matrix is too high" in str(e.value)

    with pytest.raises(ValueError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list="year", analysis_perperiod="M")
    assert "all nans" in str(e.value)
    
    with pytest.raises(ValueError) as e:
        _ = hatyan.analysis(ts=ts_pd.iloc[:1], const_list="year")
    assert "provided timeseries is less than 2 timesteps long" in str(e.value)
    
    ts_pd_manynans = ts_pd.copy()
    ts_pd_manynans.iloc[1:] = np.nan
    with pytest.raises(ValueError) as e:
        _ = hatyan.analysis(ts=ts_pd_manynans, const_list="year")
    assert "provided timeseries is less than 2 timesteps long" in str(e.value)

    ts_pd_manynans = ts_pd.copy()
    ts_pd_manynans.iloc[2:] = np.nan
    with pytest.raises(MatrixConditionTooHigh) as e:
        _ = hatyan.analysis(ts=ts_pd_manynans, const_list="year")
    assert "condition of xTx matrix is too high" in str(e.value)


@pytest.mark.unittest
def test_analysis_invalid_timeseries_index():
    file_dia = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_pd = hatyan.read_dia(file_dia)
    
    ts_pd_wrongindex = ts_pd.reset_index()
    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd_wrongindex, const_list="year")
    assert "ts.index is not of expected type" in str(e.value)
    
    ts_pd_duplicated = pd.concat([ts_pd,ts_pd])
    with pytest.raises(ValueError) as e:
        _ = hatyan.analysis(ts=ts_pd_duplicated, const_list="year")
    assert "duplicate timesteps in provided timeseries" in str(e.value)


@pytest.mark.unittest
def test_analysis_settings_invalid_values():
    file_dia = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_pd = hatyan.read_dia(file_dia)
    
    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], nodalfactors=1)
    assert str(e.value) == "invalid nodalfactors=1 type, should be bool"

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], xfac=1)
    assert str(e.value) == "invalid xfac=1 type, should be bool or dict"

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], fu_alltimes=1)
    assert str(e.value) == "invalid fu_alltimes=1 type, should be bool"

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], source=1)
    assert str(e.value) == "invalid source=1 type, should be str"

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], source="aa")
    assert str(e.value) == 'invalid source aa, should be "schureman" or "foreman"'

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], analysis_perperiod="aa")
    assert str(e.value) == 'invalid analysis_perperiod=aa type, should be False or Y/Q/M'

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], return_allperiods=1)
    assert str(e.value) == 'invalid return_allperiods=1 type, should be bool'

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], return_allperiods=True)
    assert str(e.value) == (
        'return_allperiods=True, but analysis_perperiod=False, this is not supported.'
        )

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], cs_comps=1)
    assert str(e.value) == "invalid cs_comps type, should be dict"

    with pytest.raises(KeyError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], cs_comps={})
    assert str(e.value) == "'cs_comps does not contain CS_comps_derive'"

    with pytest.raises(TypeError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2"], max_matrix_condition="aa")
    assert str(e.value) == "invalid aa type, should be int or float"


@pytest.mark.unittest
def test_analysis_invalid_constituent_lists():
    file_dia = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_pd = hatyan.read_dia(file_dia)
    
    with pytest.raises(KeyError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list="M2")
    assert 'listtype "M2" is not implemented in hatyan.get_const_list_hatyan' in str(e.value)

    with pytest.raises(ValueError) as e:
        _ = hatyan.analysis(ts=ts_pd, const_list=["M2","M2"])
    assert "remove duplicate constituents from const_list" in str(e.value)


@pytest.mark.unittest
def test_analysis_component_splitting_invalid_components():
    file_dia = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_pd = hatyan.read_dia(file_dia)
    cs_comps_valid = pd.DataFrame({'CS_comps_derive':['P1','NU2','LABDA2','K2','T2'],
                             'CS_comps_from':['K1','N2','2MN2','S2','S2'],
                             'CS_ampfacs':[0.36,0.38,0.44,0.30,0.07],
                             'CS_degincrs':[5,-22,180,3,-22]})
    with pytest.raises(ValueError) as e:
        _= hatyan.analysis(ts=ts_pd, const_list="year", cs_comps=cs_comps_valid)
    assert str(e.value) == "component P1 requested via component splitting, but already present in const_list"


@pytest.mark.unittest
def test_analysis_deprecatedsettings():
    times_pred = slice(dt.datetime(2009,12,31,14),dt.datetime(2010,1,2,12), "1min")
    file_data_comp0 = os.path.join(dir_testdata,'DENHDR_ana.txt')
    comp = hatyan.read_components(filename=file_data_comp0)
    
    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.analysis(comp=comp, times=times_pred, CS_comps=1)
    assert str(e.value) == "Argument 'CS_comps' has been deprecated for hatyan.analysis(), use 'cs_comps' instead"

    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.analysis(comp=comp, times=times_pred, xTxmat_condition_max=1)
    assert str(e.value) == "Argument 'xTxmat_condition_max' has been deprecated for hatyan.analysis(), use 'max_matrix_condition' instead"


@pytest.mark.unittest
def test_prediction_settings():
    times_pred = slice(dt.datetime(2009,12,31,14),dt.datetime(2010,1,2,12), "1min")
    current_station = 'DENHDR'
    
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    COMP_merged.attrs["nodalfactors"] = True
    COMP_merged.attrs["fu_alltimes"] = False
    COMP_merged.attrs["xfac"] = True
    ts_pred = hatyan.prediction(comp=COMP_merged, times=times_pred)
    expected = np.array([-0.55761625, -0.55022362, -0.54263955, -0.53486334, -0.52689445,
       -0.51873254, -0.51037744, -0.50182919, -0.49308803, -0.48415438])
    assert np.allclose(ts_pred["values"].values[:10], expected)
    
    COMP_merged.attrs["nodalfactors"] = True
    COMP_merged.attrs["fu_alltimes"] = True
    COMP_merged.attrs["xfac"] = True
    ts_pred = hatyan.prediction(comp=COMP_merged, times=times_pred)
    expected = np.array([-0.55758197, -0.55018925, -0.54260512, -0.53482887, -0.52685996,
           -0.51869805, -0.51034297, -0.50179476, -0.49305366, -0.4841201 ])
    assert np.allclose(ts_pred["values"].values[:10], expected)
    
    COMP_merged.attrs["nodalfactors"] = True
    COMP_merged.attrs["fu_alltimes"] = True
    COMP_merged.attrs["xfac"] = {"M2":0.8}
    ts_pred = hatyan.prediction(comp=COMP_merged, times=times_pred)
    expected = np.array([-0.55220418, -0.54482864, -0.53726402, -0.52950961, -0.52156487,
           -0.51342947, -0.50510324, -0.49658621, -0.48787861, -0.47898087])
    assert np.allclose(ts_pred["values"].values[:10], expected)
    
    COMP_merged.attrs["nodalfactors"] = True
    COMP_merged.attrs["fu_alltimes"] = True
    COMP_merged.attrs["xfac"] = False
    ts_pred = hatyan.prediction(comp=COMP_merged, times=times_pred)
    expected = np.array([-0.55115935, -0.54379567, -0.53624298, -0.52850057, -0.52056791,
           -0.51244465, -0.50413063, -0.49562589, -0.48693065, -0.47804533])
    assert np.allclose(ts_pred["values"].values[:10], expected)
    
    COMP_merged.attrs["nodalfactors"] = True
    COMP_merged.attrs["fu_alltimes"] = False
    COMP_merged.attrs["xfac"] = False
    ts_pred = hatyan.prediction(comp=COMP_merged, times=times_pred)
    expected = np.array([-0.55120816, -0.54384448, -0.53629177, -0.52854932, -0.52061659,
           -0.51249324, -0.5041791 , -0.49567421, -0.4869788 , -0.47809328])
    assert np.allclose(ts_pred["values"].values[:10], expected)
    
    COMP_merged.attrs["nodalfactors"] = False
    COMP_merged.attrs["fu_alltimes"] = False
    COMP_merged.attrs["xfac"] = False
    ts_pred = hatyan.prediction(comp=COMP_merged, times=times_pred)
    expected = np.array([-0.61931214, -0.61234662, -0.60517618, -0.59779951, -0.59021548,
           -0.58242314, -0.57442171, -0.56621065, -0.55778958, -0.54915834])
    assert np.allclose(ts_pred["values"].values[:10], expected)


@pytest.mark.unittest
def test_prediction_settings_invalid_values():
    current_station = 'VLISSGN'
    file_comp = os.path.join(dir_testdata, f'{current_station}_ana.txt')
    
    # missing xfac attribute
    comp = hatyan.read_components(filename=file_comp)
    comp.attrs.pop("xfac")
    with pytest.raises(KeyError) as e:
        hatyan.prediction(comp)
    assert str(e.value) == "'xfac'"
    
    # no settings attrs at all
    comp = hatyan.read_components(filename=file_comp)
    comp.attrs = {}
    with pytest.raises(KeyError) as e:
        hatyan.prediction(comp)
    assert str(e.value) == "'nodalfactors'"

    # incorrect settings
    comp = hatyan.read_components(filename=file_comp)
    comp.attrs["nodalfactors"] = 1
    with pytest.raises(TypeError) as e:
        hatyan.prediction(comp)
    assert str(e.value) == "invalid nodalfactors=1 type, should be bool"

    comp = hatyan.read_components(filename=file_comp)
    comp.attrs["fu_alltimes"] = 1
    with pytest.raises(TypeError) as e:
        hatyan.prediction(comp)
    assert str(e.value) == "invalid fu_alltimes=1 type, should be bool"

    comp = hatyan.read_components(filename=file_comp)
    comp.attrs["xfac"] = 1
    with pytest.raises(TypeError) as e:
        hatyan.prediction(comp)
    assert str(e.value) == "invalid xfac=1 type, should be bool or dict"

    comp = hatyan.read_components(filename=file_comp)
    comp.attrs["source"] = "aa"
    with pytest.raises(TypeError) as e:
        hatyan.prediction(comp)
    assert str(e.value) == 'invalid source aa, should be "schureman" or "foreman"'


@pytest.mark.unittest
def test_prediction_settings_invalid_times():
    current_station = 'VLISSGN'
    file_comp = os.path.join(dir_testdata, f'{current_station}_ana.txt')
    comp = hatyan.read_components(filename=file_comp)
    
    with pytest.raises(TypeError) as e:
        hatyan.prediction(comp, times=1)
    assert "times argument can be of type" in str(e.value)
    
    with pytest.raises(TypeError) as e:
        hatyan.prediction(comp, timestep=1)
    assert str(e.value) == "prediction() atonce, so 'timestep' argument not allowed"


@pytest.mark.unittest
def test_prediction_perperiod_settings_invalid_timestep():
    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0)
    
    _, comp_all = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, 
                                          analysis_perperiod="M", return_allperiods=True)
    
    with pytest.raises(ValueError) as e:
        hatyan.prediction(comp=comp_all, timestep=60)
    assert str(e.value) == "Invalid frequency: 60"

    with pytest.raises(TypeError) as e:
        hatyan.prediction(comp=comp_all, timestep=60, times=60)
    assert "prediction() per period, so 'times' argument not allowed" == str(e.value)


@pytest.mark.unittest
def test_analysis_1018():
    """
    this tests checks if we can do an analysis with timesteps that are outofbounds 
    in pandas when using unit="ns" (nanoseconds). 
    We replace the times from a dataframe with these timesteps, but with unit="us" (microseconds). 
    The analysis succeeds with pandas versions that support the unit argument, so all versions that hatyan supports.
    However, if a diafile would contain 1018 timestamps, the process will fail since 
    the dia reader will still read it as nanoseconds. This might also be the case 
    when writing the prediction to a diafile.
    """
    file_dia = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_pd = hatyan.read_dia(file_dia)
    
    # convert to year 1018
    ts_pd.index = pd.date_range('1018-01-01 00:00:00+01:00', '1018-12-31 23:00:00+01:00', freq="60min", unit="us")
    
    #prediction
    comp = hatyan.analysis(ts=ts_pd, const_list="year")
    
    assert np.allclose(comp.loc["M2"].values, [ 1.78668261, 45.85671448])


@pytest.mark.unittest
def test_prediction_1018():
    current_station = 'DENHDR'
    
    times_pred = slice(dt.datetime(1018,7,21),dt.datetime(1018,7,21,3), "10min")
    
    file_data_comp0 = os.path.join(dir_testdata,'%s_ana.txt'%(current_station))
    COMP_merged = hatyan.read_components(filename=file_data_comp0)
    
    #prediction
    ts_prediction = hatyan.prediction(comp=COMP_merged, times=times_pred)
    
    ts_prediction_times = np.array([dt.datetime(1018, 7, 21, 0, 0),
                                    dt.datetime(1018, 7, 21, 0, 10),
                                    dt.datetime(1018, 7, 21, 0, 20),
                                    dt.datetime(1018, 7, 21, 0, 30),
                                    dt.datetime(1018, 7, 21, 0, 40),
                                    dt.datetime(1018, 7, 21, 0, 50),
                                    dt.datetime(1018, 7, 21, 1, 0),
                                    dt.datetime(1018, 7, 21, 1, 10),
                                    dt.datetime(1018, 7, 21, 1, 20),
                                    dt.datetime(1018, 7, 21, 1, 30),
                                    dt.datetime(1018, 7, 21, 1, 40),
                                    dt.datetime(1018, 7, 21, 1, 50),
                                    dt.datetime(1018, 7, 21, 2, 0),
                                    dt.datetime(1018, 7, 21, 2, 10),
                                    dt.datetime(1018, 7, 21, 2, 20),
                                    dt.datetime(1018, 7, 21, 2, 30),
                                    dt.datetime(1018, 7, 21, 2, 40),
                                    dt.datetime(1018, 7, 21, 2, 50),
                                    dt.datetime(1018, 7, 21, 3, 0)], dtype=np.datetime64)
    ts_prediction_times_pd = pd.DatetimeIndex(ts_prediction_times).tz_localize("UTC+01:00")
    ts_prediction_vals = np.array([-0.77416669, -0.78821831, -0.79718123, -0.80126649, -0.80055319,
                                   -0.79466499, -0.78252257, -0.76226094, -0.73136835, -0.68705389,
                                   -0.62679467, -0.54896408, -0.4534113 , -0.34185962, -0.218017  ,
                                   -0.08734459,  0.04350389,  0.16749106,  0.27811846])
    
    assert (np.abs(ts_prediction['values'].values-ts_prediction_vals) < 10E-9).all()
    assert (np.abs(ts_prediction.index-ts_prediction_times_pd) < dt.timedelta(days=10E-9)).all()


@pytest.mark.unittest
def test_prediction_comp_and_times_different_timezones():
    """
    it is possible to supply a timezone via times, but the comp dataframe also contains a timezone already.
    The components timezone is leading, but the times will be converted to that timezone also.
    From the below test, both predictions are therefore UTC+01:00, but the times are shifted
    https://github.com/Deltares/hatyan/issues/334
    """
    
    current_station = 'VLISSGN'
    
    const_list = hatyan.get_const_list_hatyan('year') # 94 constituents
    
    file_comp = os.path.join(dir_testdata,f'{current_station}_obs1.txt')
    ts_meas = hatyan.read_dia(filename=file_comp, station=current_station)
    
    comp_naive = hatyan.analysis(ts=ts_meas, const_list=const_list)
    comp_naive.attrs["tzone"] = None
    comp_met = hatyan.analysis(ts=ts_meas, const_list=const_list)
    
    times_naive = slice("2019-01-01","2020-01-01", "10min")
    times_met = pd.date_range("2019-01-01","2020-01-01", freq="10min", tz="UTC+01:00")
    times_utc = pd.date_range("2019-01-01","2020-01-01", freq="10min", tz="UTC+00:00")
    
    pred_naive = hatyan.prediction(comp=comp_naive, times=times_naive)
    pred_met = hatyan.prediction(comp=comp_met, times=times_met)
    pred_utc = hatyan.prediction(comp=comp_met, times=times_utc)
    
    assert pred_naive.index.tz is None
    assert pred_naive.index.freq is not None
    assert pred_naive.index[0] == pd.Timestamp('2019-01-01 00:00:00')
    assert pred_met.index.tz == dt.timezone(dt.timedelta(seconds=3600))
    assert pred_met.index.freq is not None
    assert pred_met.index[0] == pd.Timestamp('2019-01-01 00:00:00+0100')
    assert pred_utc.index.tz == dt.timezone.utc
    assert pred_utc.index.freq is not None
    assert pred_utc.index[0] == pd.Timestamp('2019-01-01 00:00:00+0000')
    assert ((pred_naive - pred_met.tz_localize(None)).dropna()["values"] < 1e-9).all()
    assert ((pred_utc - pred_met).dropna()["values"] < 1e-9).all()


@pytest.mark.unittest
def test_prediction_times_tznaive_comp_tzaware(caplog):
    """
    https://github.com/Deltares/hatyan/issues/334
    """
    comp = pd.DataFrame({"A": [1, 0.5, 0.2],
                         "phi_deg": [10,15,20]}, 
                        index=["M2","M4","S2"])
    comp.attrs["nodalfactors"] = True
    comp.attrs["fu_alltimes"] = True
    comp.attrs["xfac"] = False
    comp.attrs["source"] = "schureman"
    comp.attrs["tzone"] = "UTC+01:00"
    dtindex = pd.date_range("2020-01-01","2020-01-02", freq="10min")
    hatyan.prediction(comp, times=dtindex)
    warning_text = ("provided times are timezone-naive and provided components are "
                    "timezone-aware. The times are being interpreted as if they would "
                    "have the same timezone as the components: UTC+01:00")
    assert warning_text in caplog.text


@pytest.mark.unittest
def test_prediction_raise_mixed_tznaive_tzaware():
    """
    https://github.com/Deltares/hatyan/issues/334
    """
    comp = pd.DataFrame({"A": [1, 0.5, 0.2],
                         "phi_deg": [10,15,20]}, 
                        index=["M2","M4","S2"])
    comp.attrs["nodalfactors"] = True
    comp.attrs["fu_alltimes"] = True
    comp.attrs["xfac"] = False
    comp.attrs["source"] = "schureman"
    comp.attrs["tzone"] = None
    dtindex = pd.date_range("2020-01-01 00:00 +00:00","2020-01-02 00:00 +00:00", freq="10min")
    with pytest.raises(ValueError) as e:
        hatyan.prediction(comp, times=dtindex)
    assert "provided times and components should both be timezone-aware or timezone-naive, not mixed." in str(e.value)


@pytest.mark.unittest
def test_prediction_comp_no_timezone():
    """
    https://github.com/Deltares/hatyan/issues/317
    """
    comp = pd.DataFrame({"A": [1, 0.5, 0.2],
                         "phi_deg": [10,15,20]}, 
                        index=["M2","M4","S2"])
    comp.attrs["nodalfactors"] = True
    comp.attrs["fu_alltimes"] = True
    comp.attrs["xfac"] = False
    comp.attrs["source"] = "schureman"
    dtindex = pd.date_range("2020-01-01","2020-01-02", freq="10min")
    pred = hatyan.prediction(comp, times=dtindex)
    assert pred.index.tz is None
    assert pred.index[0] == pd.Timestamp('2020-01-01 00:00:00')
    assert pred.index[-1] == pd.Timestamp('2020-01-02 00:00:00')


@pytest.mark.unittest
def test_prediction_perperiod_month():
    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0)
    
    comp_mean, comp_all = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, 
                                          analysis_perperiod="M", return_allperiods=True)
    
    ts_pred_atonce = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index)
    
    ts_pred_allmonths = hatyan.prediction(comp=comp_all, timestep="60min")
    
    expected_atonce = np.array([-1.39664945, -0.85846188, -0.22419191,  0.7596607 ,  1.91604743,
            2.17294986,  1.57508019,  0.88085591, -0.01977715, -1.01827975])
    expected_allmonths = np.array([-1.51643013, -1.06743931, -0.51656439,  0.40987597,  1.66862452,
            2.09709427,  1.60408414,  0.9891484 ,  0.13523824, -0.89396145])
    
    assert np.allclose(ts_pred_atonce["values"].values[:10], expected_atonce)
    assert np.allclose(ts_pred_allmonths["values"].values[:10], expected_allmonths)
    assert len(ts_pred_atonce) == len(ts_pred_allmonths)
    assert ts_pred_allmonths.index.tz == pytz.FixedOffset(60)


@pytest.mark.unittest
def test_prediction_perperiod_year():
    # TODO: when the prediction function does not treat Y and M differently, this testcase can be dropped
    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0)
    
    comp_mean, comp_all = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, 
                                          analysis_perperiod="Y", return_allperiods=True)
    
    ts_pred_atonce = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index)
    
    ts_pred_allyears = hatyan.prediction(comp=comp_all, timestep="60min")
    
    expected_atonce = np.array([-1.38346859, -0.84555223, -0.21307598,  0.76601685,  1.91438215,
            2.16373654,  1.56241672,  0.86686561, -0.0329582 , -1.02668776])
    expected_allyears = np.array([-1.38346859, -0.84555223, -0.21307598,  0.76601685,  1.91438215,
            2.16373654,  1.56241672,  0.86686561, -0.0329582 , -1.02668776])
    
    assert np.allclose(ts_pred_atonce["values"].values[:10], expected_atonce)
    assert np.allclose(ts_pred_allyears["values"].values[:10], expected_allyears)
    assert len(ts_pred_atonce) == len(ts_pred_allyears)
    assert ts_pred_allyears.index.tz == pytz.FixedOffset(60)


@pytest.mark.unittest
def test_prediction_hasrequiredargs():
    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0)
    
    comp_mean, comp_all = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True, 
                                          analysis_perperiod="M", return_allperiods=True)
    
    with pytest.raises(TypeError) as e:
        _ = hatyan.prediction(comp=comp_mean)
    assert str(e.value) == "prediction() atonce, so 'times' argument should not be None"
    
    with pytest.raises(TypeError) as e:
        _ = hatyan.prediction(comp=comp_all)
    assert str(e.value) == "prediction() per period, so 'timestep' argument should not be None"


@pytest.mark.unittest
def test_prediction_deprecatedsettings():
    file_data_comp0 = os.path.join(dir_testdata,'VLISSGN_obs1.txt')
    ts_measurements_group0 = hatyan.read_dia(filename=file_data_comp0)
    
    comp_mean = hatyan.analysis(ts=ts_measurements_group0, const_list='month', nodalfactors=True, fu_alltimes=False, xfac=True)
    
    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index, times_pred_all="")
    assert str(e.value) == "Argument 'times_pred_all' has been deprecated for hatyan.prediction(), use 'times' instead"
    
    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index, times_ext="")
    assert str(e.value) == "Argument 'times_ext' has been deprecated for hatyan.prediction(), use 'times' instead"
    
    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index, nodalfactors=False)
    assert "as attribute of the component dataframe" in str(e.value)
        
    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index, xfac=False)
    assert "as attribute of the component dataframe" in str(e.value)
    
    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index, fu_alltimes=False)
    assert "as attribute of the component dataframe" in str(e.value)
    
    with pytest.raises(DeprecationWarning) as e:
        _ = hatyan.prediction(comp=comp_mean, times=ts_measurements_group0.index, source=False)
    assert "as attribute of the component dataframe" in str(e.value)
