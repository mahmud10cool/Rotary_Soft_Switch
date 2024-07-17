    function targMap = targDataMap(),

    ;%***********************
    ;% Create Parameter Map *
    ;%***********************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 3;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc paramMap
        ;%
        paramMap.nSections           = nTotSects;
        paramMap.sectIdxOffset       = sectIdxOffset;
            paramMap.sections(nTotSects) = dumSection; %prealloc
        paramMap.nTotData            = -1;

        ;%
        ;% Auto data (rtP)
        ;%
            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% rtP.param
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            paramMap.sections(1) = section;
            clear section

            section.nData     = 44;
            section.data(44)  = dumData; %prealloc

                    ;% rtP.Out1_Y0
                    section.data(1).logicalSrcIdx = 1;
                    section.data(1).dtTransOffset = 0;

                    ;% rtP.Integrator4_IC
                    section.data(2).logicalSrcIdx = 2;
                    section.data(2).dtTransOffset = 1;

                    ;% rtP.Gain_Gain
                    section.data(3).logicalSrcIdx = 3;
                    section.data(3).dtTransOffset = 2;

                    ;% rtP.Relay_OnVal
                    section.data(4).logicalSrcIdx = 4;
                    section.data(4).dtTransOffset = 3;

                    ;% rtP.Relay_OffVal
                    section.data(5).logicalSrcIdx = 5;
                    section.data(5).dtTransOffset = 4;

                    ;% rtP.Relay_YOn
                    section.data(6).logicalSrcIdx = 6;
                    section.data(6).dtTransOffset = 5;

                    ;% rtP.Relay_YOff
                    section.data(7).logicalSrcIdx = 7;
                    section.data(7).dtTransOffset = 6;

                    ;% rtP.Integrator13_IC
                    section.data(8).logicalSrcIdx = 8;
                    section.data(8).dtTransOffset = 7;

                    ;% rtP.Integrator14_IC
                    section.data(9).logicalSrcIdx = 9;
                    section.data(9).dtTransOffset = 8;

                    ;% rtP.Integrator7_IC
                    section.data(10).logicalSrcIdx = 10;
                    section.data(10).dtTransOffset = 9;

                    ;% rtP.Integrator12_IC
                    section.data(11).logicalSrcIdx = 11;
                    section.data(11).dtTransOffset = 10;

                    ;% rtP.Integrator10_IC
                    section.data(12).logicalSrcIdx = 12;
                    section.data(12).dtTransOffset = 11;

                    ;% rtP.Integrator8_IC
                    section.data(13).logicalSrcIdx = 13;
                    section.data(13).dtTransOffset = 12;

                    ;% rtP.Integrator9_IC
                    section.data(14).logicalSrcIdx = 14;
                    section.data(14).dtTransOffset = 13;

                    ;% rtP.Integrator11_IC
                    section.data(15).logicalSrcIdx = 15;
                    section.data(15).dtTransOffset = 14;

                    ;% rtP.Gain3_Gain
                    section.data(16).logicalSrcIdx = 16;
                    section.data(16).dtTransOffset = 15;

                    ;% rtP.Integrator6_IC
                    section.data(17).logicalSrcIdx = 17;
                    section.data(17).dtTransOffset = 16;

                    ;% rtP.Integrator16_IC
                    section.data(18).logicalSrcIdx = 18;
                    section.data(18).dtTransOffset = 17;

                    ;% rtP.TransferFcn_A
                    section.data(19).logicalSrcIdx = 19;
                    section.data(19).dtTransOffset = 18;

                    ;% rtP.TransferFcn_C
                    section.data(20).logicalSrcIdx = 20;
                    section.data(20).dtTransOffset = 20;

                    ;% rtP.Saturation_LowerSat
                    section.data(21).logicalSrcIdx = 21;
                    section.data(21).dtTransOffset = 22;

                    ;% rtP.TransferFcn_A_msxihde20w
                    section.data(22).logicalSrcIdx = 22;
                    section.data(22).dtTransOffset = 23;

                    ;% rtP.TransferFcn_C_nnurr522qa
                    section.data(23).logicalSrcIdx = 23;
                    section.data(23).dtTransOffset = 25;

                    ;% rtP.Saturation_LowerSat_l0ksmwxpmq
                    section.data(24).logicalSrcIdx = 24;
                    section.data(24).dtTransOffset = 27;

                    ;% rtP.TransferFcn_A_bwuzoqq3iu
                    section.data(25).logicalSrcIdx = 25;
                    section.data(25).dtTransOffset = 28;

                    ;% rtP.TransferFcn_C_j01bq5y01n
                    section.data(26).logicalSrcIdx = 26;
                    section.data(26).dtTransOffset = 30;

                    ;% rtP.Saturation_LowerSat_m0gx2ayxie
                    section.data(27).logicalSrcIdx = 27;
                    section.data(27).dtTransOffset = 32;

                    ;% rtP.TransferFcn_A_cqt2tjfue3
                    section.data(28).logicalSrcIdx = 28;
                    section.data(28).dtTransOffset = 33;

                    ;% rtP.TransferFcn_C_kavhf31tp1
                    section.data(29).logicalSrcIdx = 29;
                    section.data(29).dtTransOffset = 35;

                    ;% rtP.Saturation_LowerSat_l5sm1vdo10
                    section.data(30).logicalSrcIdx = 30;
                    section.data(30).dtTransOffset = 37;

                    ;% rtP.Integrator5_IC
                    section.data(31).logicalSrcIdx = 31;
                    section.data(31).dtTransOffset = 38;

                    ;% rtP.Delay_InitialCondition
                    section.data(32).logicalSrcIdx = 32;
                    section.data(32).dtTransOffset = 39;

                    ;% rtP.Velocity_Amp
                    section.data(33).logicalSrcIdx = 33;
                    section.data(33).dtTransOffset = 40;

                    ;% rtP.Velocity_Bias
                    section.data(34).logicalSrcIdx = 34;
                    section.data(34).dtTransOffset = 41;

                    ;% rtP.Velocity_Freq
                    section.data(35).logicalSrcIdx = 35;
                    section.data(35).dtTransOffset = 42;

                    ;% rtP.Velocity_Phase
                    section.data(36).logicalSrcIdx = 36;
                    section.data(36).dtTransOffset = 43;

                    ;% rtP.Integrator15_IC
                    section.data(37).logicalSrcIdx = 37;
                    section.data(37).dtTransOffset = 44;

                    ;% rtP.TransferFcn_A_m2qs34lzon
                    section.data(38).logicalSrcIdx = 38;
                    section.data(38).dtTransOffset = 45;

                    ;% rtP.TransferFcn_C_lu1wzht1on
                    section.data(39).logicalSrcIdx = 39;
                    section.data(39).dtTransOffset = 47;

                    ;% rtP.Saturation_LowerSat_mclkng5niv
                    section.data(40).logicalSrcIdx = 40;
                    section.data(40).dtTransOffset = 49;

                    ;% rtP.InitialSpeed_Value
                    section.data(41).logicalSrcIdx = 41;
                    section.data(41).dtTransOffset = 50;

                    ;% rtP.Gain5_Gain
                    section.data(42).logicalSrcIdx = 42;
                    section.data(42).dtTransOffset = 51;

                    ;% rtP.Constant_Value
                    section.data(43).logicalSrcIdx = 43;
                    section.data(43).dtTransOffset = 52;

                    ;% rtP.Constant5_Value
                    section.data(44).logicalSrcIdx = 44;
                    section.data(44).dtTransOffset = 53;

            nTotData = nTotData + section.nData;
            paramMap.sections(2) = section;
            clear section

            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% rtP.ManualSwitch_CurrentSetting
                    section.data(1).logicalSrcIdx = 45;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            paramMap.sections(3) = section;
            clear section


            ;%
            ;% Non-auto Data (parameter)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        paramMap.nTotData = nTotData;



    ;%**************************
    ;% Create Block Output Map *
    ;%**************************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 1;
        sectIdxOffset = 0;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc sigMap
        ;%
        sigMap.nSections           = nTotSects;
        sigMap.sectIdxOffset       = sectIdxOffset;
            sigMap.sections(nTotSects) = dumSection; %prealloc
        sigMap.nTotData            = -1;

        ;%
        ;% Auto data (rtB)
        ;%
            section.nData     = 79;
            section.data(79)  = dumData; %prealloc

                    ;% rtB.du3hz2fyzs
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% rtB.frh4wmn3vf
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

                    ;% rtB.mcwb4pv1pd
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 2;

                    ;% rtB.cqrrgwzxyv
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 3;

                    ;% rtB.dxeimajz1j
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 4;

                    ;% rtB.a5ddnr2k0h
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 5;

                    ;% rtB.blk35wfutl
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 6;

                    ;% rtB.ikogzwrhkz
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 7;

                    ;% rtB.fakrdhsmg5
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 8;

                    ;% rtB.cy5jq00nod
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 9;

                    ;% rtB.ht3vr3fcw1
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 10;

                    ;% rtB.axh1siq1cq
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 11;

                    ;% rtB.fkszch0qwe
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 12;

                    ;% rtB.idqoyikyv2
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 13;

                    ;% rtB.hgyc45hfyg
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 14;

                    ;% rtB.ku0vq41nvc
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 15;

                    ;% rtB.cjeggiimao
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 16;

                    ;% rtB.ab5hqtcfdy
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 17;

                    ;% rtB.ip00peawg3
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 18;

                    ;% rtB.iubu1mncgz
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 19;

                    ;% rtB.na1kfpe1mt
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 20;

                    ;% rtB.eee41maadj
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 21;

                    ;% rtB.pebffqxdnd
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 22;

                    ;% rtB.p2rr4uzpkw
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 23;

                    ;% rtB.fj3bxcs4ue
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 24;

                    ;% rtB.kkbkkd5z5m
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 25;

                    ;% rtB.gztmpcdsdg
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 26;

                    ;% rtB.mxnffxyshc
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 27;

                    ;% rtB.fr14joswus
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 28;

                    ;% rtB.g3vckbtrp3
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 29;

                    ;% rtB.cqkdaw1iu3
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 30;

                    ;% rtB.o123mdetox
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 31;

                    ;% rtB.jsgd1kggkx
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 32;

                    ;% rtB.nc3ujqp1pg
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 33;

                    ;% rtB.iilt4zdhq4
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 34;

                    ;% rtB.foabxlso54
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 35;

                    ;% rtB.mscp3mtqiu
                    section.data(37).logicalSrcIdx = 36;
                    section.data(37).dtTransOffset = 36;

                    ;% rtB.orhxqsmvfm
                    section.data(38).logicalSrcIdx = 37;
                    section.data(38).dtTransOffset = 37;

                    ;% rtB.pcceognrbb
                    section.data(39).logicalSrcIdx = 38;
                    section.data(39).dtTransOffset = 38;

                    ;% rtB.hcfgbzdbcf
                    section.data(40).logicalSrcIdx = 39;
                    section.data(40).dtTransOffset = 39;

                    ;% rtB.jost0dydiu
                    section.data(41).logicalSrcIdx = 40;
                    section.data(41).dtTransOffset = 40;

                    ;% rtB.nodqar1hqn
                    section.data(42).logicalSrcIdx = 41;
                    section.data(42).dtTransOffset = 41;

                    ;% rtB.oi30omicce
                    section.data(43).logicalSrcIdx = 42;
                    section.data(43).dtTransOffset = 42;

                    ;% rtB.nmkcufy13s
                    section.data(44).logicalSrcIdx = 43;
                    section.data(44).dtTransOffset = 43;

                    ;% rtB.bhydlzvdgu
                    section.data(45).logicalSrcIdx = 44;
                    section.data(45).dtTransOffset = 44;

                    ;% rtB.jnk41mcgi4
                    section.data(46).logicalSrcIdx = 45;
                    section.data(46).dtTransOffset = 45;

                    ;% rtB.nhtjb14mqm
                    section.data(47).logicalSrcIdx = 46;
                    section.data(47).dtTransOffset = 46;

                    ;% rtB.becosfp3hi
                    section.data(48).logicalSrcIdx = 47;
                    section.data(48).dtTransOffset = 47;

                    ;% rtB.lphlstzov0
                    section.data(49).logicalSrcIdx = 48;
                    section.data(49).dtTransOffset = 48;

                    ;% rtB.dcpeged5b2
                    section.data(50).logicalSrcIdx = 49;
                    section.data(50).dtTransOffset = 49;

                    ;% rtB.eqjfsu4kur
                    section.data(51).logicalSrcIdx = 50;
                    section.data(51).dtTransOffset = 50;

                    ;% rtB.cpcl24qdbb
                    section.data(52).logicalSrcIdx = 51;
                    section.data(52).dtTransOffset = 51;

                    ;% rtB.hpzll4vow0
                    section.data(53).logicalSrcIdx = 52;
                    section.data(53).dtTransOffset = 52;

                    ;% rtB.okpqgwd3h2
                    section.data(54).logicalSrcIdx = 53;
                    section.data(54).dtTransOffset = 53;

                    ;% rtB.cbpdlgxiw5
                    section.data(55).logicalSrcIdx = 54;
                    section.data(55).dtTransOffset = 54;

                    ;% rtB.lv3onhl4x1
                    section.data(56).logicalSrcIdx = 55;
                    section.data(56).dtTransOffset = 55;

                    ;% rtB.fmn4xw0byt
                    section.data(57).logicalSrcIdx = 56;
                    section.data(57).dtTransOffset = 56;

                    ;% rtB.aeifxkbw2r
                    section.data(58).logicalSrcIdx = 57;
                    section.data(58).dtTransOffset = 57;

                    ;% rtB.idpdlwl5kt
                    section.data(59).logicalSrcIdx = 58;
                    section.data(59).dtTransOffset = 58;

                    ;% rtB.haqtqgqlej
                    section.data(60).logicalSrcIdx = 59;
                    section.data(60).dtTransOffset = 59;

                    ;% rtB.e0pcr2qr33
                    section.data(61).logicalSrcIdx = 60;
                    section.data(61).dtTransOffset = 60;

                    ;% rtB.fo05goe5wa
                    section.data(62).logicalSrcIdx = 61;
                    section.data(62).dtTransOffset = 61;

                    ;% rtB.a1drd1cazf
                    section.data(63).logicalSrcIdx = 62;
                    section.data(63).dtTransOffset = 62;

                    ;% rtB.lqpliiebzv
                    section.data(64).logicalSrcIdx = 63;
                    section.data(64).dtTransOffset = 63;

                    ;% rtB.hrpirqqfri
                    section.data(65).logicalSrcIdx = 64;
                    section.data(65).dtTransOffset = 64;

                    ;% rtB.ollqjfiay1
                    section.data(66).logicalSrcIdx = 65;
                    section.data(66).dtTransOffset = 65;

                    ;% rtB.fzo2unh14f
                    section.data(67).logicalSrcIdx = 66;
                    section.data(67).dtTransOffset = 66;

                    ;% rtB.msv1wjj5cz
                    section.data(68).logicalSrcIdx = 67;
                    section.data(68).dtTransOffset = 67;

                    ;% rtB.lxsvfwftlj
                    section.data(69).logicalSrcIdx = 68;
                    section.data(69).dtTransOffset = 68;

                    ;% rtB.jjwnm1umsr
                    section.data(70).logicalSrcIdx = 69;
                    section.data(70).dtTransOffset = 69;

                    ;% rtB.jbbgrda5ov
                    section.data(71).logicalSrcIdx = 70;
                    section.data(71).dtTransOffset = 70;

                    ;% rtB.bojkujnbsq
                    section.data(72).logicalSrcIdx = 71;
                    section.data(72).dtTransOffset = 71;

                    ;% rtB.bnvm1inada
                    section.data(73).logicalSrcIdx = 72;
                    section.data(73).dtTransOffset = 72;

                    ;% rtB.dw11js5duv
                    section.data(74).logicalSrcIdx = 73;
                    section.data(74).dtTransOffset = 73;

                    ;% rtB.jabohzh110
                    section.data(75).logicalSrcIdx = 74;
                    section.data(75).dtTransOffset = 74;

                    ;% rtB.iikt2yfv4y
                    section.data(76).logicalSrcIdx = 75;
                    section.data(76).dtTransOffset = 75;

                    ;% rtB.ilwj52njki
                    section.data(77).logicalSrcIdx = 76;
                    section.data(77).dtTransOffset = 76;

                    ;% rtB.pwpjs2gjry
                    section.data(78).logicalSrcIdx = 77;
                    section.data(78).dtTransOffset = 77;

                    ;% rtB.b3cdqkxcij
                    section.data(79).logicalSrcIdx = 78;
                    section.data(79).dtTransOffset = 78;

            nTotData = nTotData + section.nData;
            sigMap.sections(1) = section;
            clear section


            ;%
            ;% Non-auto Data (signal)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        sigMap.nTotData = nTotData;



    ;%*******************
    ;% Create DWork Map *
    ;%*******************
    
        nTotData      = 0; %add to this count as we go
        nTotSects     = 7;
        sectIdxOffset = 1;

        ;%
        ;% Define dummy sections & preallocate arrays
        ;%
        dumSection.nData = -1;
        dumSection.data  = [];

        dumData.logicalSrcIdx = -1;
        dumData.dtTransOffset = -1;

        ;%
        ;% Init/prealloc dworkMap
        ;%
        dworkMap.nSections           = nTotSects;
        dworkMap.sectIdxOffset       = sectIdxOffset;
            dworkMap.sections(nTotSects) = dumSection; %prealloc
        dworkMap.nTotData            = -1;

        ;%
        ;% Auto data (rtDW)
        ;%
            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% rtDW.mog1ypd0na
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            dworkMap.sections(1) = section;
            clear section

            section.nData     = 20;
            section.data(20)  = dumData; %prealloc

                    ;% rtDW.mjjnlupcqc.LoggedData
                    section.data(1).logicalSrcIdx = 1;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.m5bl3o4gs5.LoggedData
                    section.data(2).logicalSrcIdx = 2;
                    section.data(2).dtTransOffset = 2;

                    ;% rtDW.pojnsaofwq.LoggedData
                    section.data(3).logicalSrcIdx = 3;
                    section.data(3).dtTransOffset = 3;

                    ;% rtDW.ajsjglxjta.LoggedData
                    section.data(4).logicalSrcIdx = 4;
                    section.data(4).dtTransOffset = 7;

                    ;% rtDW.ahihqvv0qg.LoggedData
                    section.data(5).logicalSrcIdx = 5;
                    section.data(5).dtTransOffset = 8;

                    ;% rtDW.lftwng3afb.LoggedData
                    section.data(6).logicalSrcIdx = 6;
                    section.data(6).dtTransOffset = 11;

                    ;% rtDW.dgbo00d2tw.LoggedData
                    section.data(7).logicalSrcIdx = 7;
                    section.data(7).dtTransOffset = 12;

                    ;% rtDW.cpa2ih4ymh.AQHandles
                    section.data(8).logicalSrcIdx = 8;
                    section.data(8).dtTransOffset = 13;

                    ;% rtDW.igylchmu0t.AQHandles
                    section.data(9).logicalSrcIdx = 9;
                    section.data(9).dtTransOffset = 14;

                    ;% rtDW.aqjm4jt1o5.AQHandles
                    section.data(10).logicalSrcIdx = 10;
                    section.data(10).dtTransOffset = 15;

                    ;% rtDW.gzhbrbv4c4.AQHandles
                    section.data(11).logicalSrcIdx = 11;
                    section.data(11).dtTransOffset = 16;

                    ;% rtDW.htile1lzbu.AQHandles
                    section.data(12).logicalSrcIdx = 12;
                    section.data(12).dtTransOffset = 17;

                    ;% rtDW.p1ayluai4b.AQHandles
                    section.data(13).logicalSrcIdx = 13;
                    section.data(13).dtTransOffset = 18;

                    ;% rtDW.gjtuyi3wgw.AQHandles
                    section.data(14).logicalSrcIdx = 14;
                    section.data(14).dtTransOffset = 19;

                    ;% rtDW.jd4fubbir0.AQHandles
                    section.data(15).logicalSrcIdx = 15;
                    section.data(15).dtTransOffset = 20;

                    ;% rtDW.kbuynfluan.AQHandles
                    section.data(16).logicalSrcIdx = 16;
                    section.data(16).dtTransOffset = 21;

                    ;% rtDW.j5pyh2o05c.AQHandles
                    section.data(17).logicalSrcIdx = 17;
                    section.data(17).dtTransOffset = 22;

                    ;% rtDW.p0fxbp34np.AQHandles
                    section.data(18).logicalSrcIdx = 18;
                    section.data(18).dtTransOffset = 23;

                    ;% rtDW.n4qxc3kzkl.AQHandles
                    section.data(19).logicalSrcIdx = 19;
                    section.data(19).dtTransOffset = 24;

                    ;% rtDW.mma0ea0nbs.LoggedData
                    section.data(20).logicalSrcIdx = 20;
                    section.data(20).dtTransOffset = 25;

            nTotData = nTotData + section.nData;
            dworkMap.sections(2) = section;
            clear section

            section.nData     = 16;
            section.data(16)  = dumData; %prealloc

                    ;% rtDW.d1ijg3op2k
                    section.data(1).logicalSrcIdx = 21;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.brka54cpul
                    section.data(2).logicalSrcIdx = 22;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.ojcciswcd3
                    section.data(3).logicalSrcIdx = 23;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.nfkviszx2j
                    section.data(4).logicalSrcIdx = 24;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.ndyrr3tae0
                    section.data(5).logicalSrcIdx = 25;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.c2jtsctdjh
                    section.data(6).logicalSrcIdx = 26;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.mevs4smziy
                    section.data(7).logicalSrcIdx = 27;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.egiikq3opb
                    section.data(8).logicalSrcIdx = 28;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.puivctwyq5
                    section.data(9).logicalSrcIdx = 29;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.bx2j5ck3fy
                    section.data(10).logicalSrcIdx = 30;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.fzispyfeqr
                    section.data(11).logicalSrcIdx = 31;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.brdpyuefrq
                    section.data(12).logicalSrcIdx = 32;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.esun0jx3xw
                    section.data(13).logicalSrcIdx = 33;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.ncj5ml0ihc
                    section.data(14).logicalSrcIdx = 34;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.kjn3oyfrac
                    section.data(15).logicalSrcIdx = 35;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.d3q0cyta13
                    section.data(16).logicalSrcIdx = 36;
                    section.data(16).dtTransOffset = 15;

            nTotData = nTotData + section.nData;
            dworkMap.sections(3) = section;
            clear section

            section.nData     = 10;
            section.data(10)  = dumData; %prealloc

                    ;% rtDW.di0xkhax5v
                    section.data(1).logicalSrcIdx = 37;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.dol3avfpna
                    section.data(2).logicalSrcIdx = 38;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.mfw5n1ea1y
                    section.data(3).logicalSrcIdx = 39;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.oj5q0jsluk
                    section.data(4).logicalSrcIdx = 40;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.dfi4io44l1
                    section.data(5).logicalSrcIdx = 41;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.bapvsmhne1
                    section.data(6).logicalSrcIdx = 42;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.gmzfgv4p33
                    section.data(7).logicalSrcIdx = 43;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.bxfqhgcloy
                    section.data(8).logicalSrcIdx = 44;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.gwwzsdaztz
                    section.data(9).logicalSrcIdx = 45;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.agicmdjqdi
                    section.data(10).logicalSrcIdx = 46;
                    section.data(10).dtTransOffset = 9;

            nTotData = nTotData + section.nData;
            dworkMap.sections(4) = section;
            clear section

            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% rtDW.astl0qyquk
                    section.data(1).logicalSrcIdx = 47;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            dworkMap.sections(5) = section;
            clear section

            section.nData     = 16;
            section.data(16)  = dumData; %prealloc

                    ;% rtDW.ekvnfwjtsb
                    section.data(1).logicalSrcIdx = 48;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.dql53suw4t
                    section.data(2).logicalSrcIdx = 49;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.gz0frs1wku
                    section.data(3).logicalSrcIdx = 50;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.iwrzpan53u
                    section.data(4).logicalSrcIdx = 51;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.ivbuycouap
                    section.data(5).logicalSrcIdx = 52;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.oatrzraqor
                    section.data(6).logicalSrcIdx = 53;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.h3yqjmxcau
                    section.data(7).logicalSrcIdx = 54;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.c2gp1hcawr
                    section.data(8).logicalSrcIdx = 55;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.kzezuc5kmj
                    section.data(9).logicalSrcIdx = 56;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.baej0t4mz4
                    section.data(10).logicalSrcIdx = 57;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.e4o1f3himj
                    section.data(11).logicalSrcIdx = 58;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.etn0cabray
                    section.data(12).logicalSrcIdx = 59;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.nq5efgjr4u
                    section.data(13).logicalSrcIdx = 60;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.ldwj0svuv5
                    section.data(14).logicalSrcIdx = 61;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.et5jt4roea
                    section.data(15).logicalSrcIdx = 62;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.jhgzfejyh1
                    section.data(16).logicalSrcIdx = 63;
                    section.data(16).dtTransOffset = 15;

            nTotData = nTotData + section.nData;
            dworkMap.sections(6) = section;
            clear section

            section.nData     = 17;
            section.data(17)  = dumData; %prealloc

                    ;% rtDW.e3b1amt0gc
                    section.data(1).logicalSrcIdx = 64;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.g3swjv3kh4
                    section.data(2).logicalSrcIdx = 65;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.cifi3jtean
                    section.data(3).logicalSrcIdx = 66;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.mx1mzmmudm
                    section.data(4).logicalSrcIdx = 67;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.nhmw0vweas
                    section.data(5).logicalSrcIdx = 68;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.enb1yfnvht
                    section.data(6).logicalSrcIdx = 69;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.au1n1jxucn
                    section.data(7).logicalSrcIdx = 70;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.cwu3pvv3pv
                    section.data(8).logicalSrcIdx = 71;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.dwzumrr4bn
                    section.data(9).logicalSrcIdx = 72;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.jelszbraal
                    section.data(10).logicalSrcIdx = 73;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.cbdxrybiau
                    section.data(11).logicalSrcIdx = 74;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.dia2paiwqh
                    section.data(12).logicalSrcIdx = 75;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.f40yghkee4
                    section.data(13).logicalSrcIdx = 76;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.pke0ukp1xy
                    section.data(14).logicalSrcIdx = 77;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.ob1pklob1o
                    section.data(15).logicalSrcIdx = 78;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.jvf3pz3ask
                    section.data(16).logicalSrcIdx = 79;
                    section.data(16).dtTransOffset = 15;

                    ;% rtDW.cy4ormxj51
                    section.data(17).logicalSrcIdx = 80;
                    section.data(17).dtTransOffset = 16;

            nTotData = nTotData + section.nData;
            dworkMap.sections(7) = section;
            clear section


            ;%
            ;% Non-auto Data (dwork)
            ;%


        ;%
        ;% Add final counts to struct.
        ;%
        dworkMap.nTotData = nTotData;



    ;%
    ;% Add individual maps to base struct.
    ;%

    targMap.paramMap  = paramMap;
    targMap.signalMap = sigMap;
    targMap.dworkMap  = dworkMap;

    ;%
    ;% Add checksums to base struct.
    ;%


    targMap.checksum0 = 2806597705;
    targMap.checksum1 = 165651654;
    targMap.checksum2 = 2791374129;
    targMap.checksum3 = 1143219654;

