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

                    ;% rtB.mzfn3lhia4
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

                    ;% rtB.jswgtkntlr
                    section.data(2).logicalSrcIdx = 1;
                    section.data(2).dtTransOffset = 1;

                    ;% rtB.pyky0iaydb
                    section.data(3).logicalSrcIdx = 2;
                    section.data(3).dtTransOffset = 2;

                    ;% rtB.cwl31aims5
                    section.data(4).logicalSrcIdx = 3;
                    section.data(4).dtTransOffset = 3;

                    ;% rtB.dk50fnzsks
                    section.data(5).logicalSrcIdx = 4;
                    section.data(5).dtTransOffset = 4;

                    ;% rtB.lr3evndp05
                    section.data(6).logicalSrcIdx = 5;
                    section.data(6).dtTransOffset = 5;

                    ;% rtB.gaymsnblvn
                    section.data(7).logicalSrcIdx = 6;
                    section.data(7).dtTransOffset = 6;

                    ;% rtB.hfumqnwosi
                    section.data(8).logicalSrcIdx = 7;
                    section.data(8).dtTransOffset = 7;

                    ;% rtB.jcnaw1aobb
                    section.data(9).logicalSrcIdx = 8;
                    section.data(9).dtTransOffset = 8;

                    ;% rtB.c5jlsr0y2m
                    section.data(10).logicalSrcIdx = 9;
                    section.data(10).dtTransOffset = 9;

                    ;% rtB.gjajnpeavu
                    section.data(11).logicalSrcIdx = 10;
                    section.data(11).dtTransOffset = 10;

                    ;% rtB.onszijeutj
                    section.data(12).logicalSrcIdx = 11;
                    section.data(12).dtTransOffset = 11;

                    ;% rtB.o5jgaxymjd
                    section.data(13).logicalSrcIdx = 12;
                    section.data(13).dtTransOffset = 12;

                    ;% rtB.aqvwe35511
                    section.data(14).logicalSrcIdx = 13;
                    section.data(14).dtTransOffset = 13;

                    ;% rtB.ga4spwnzmb
                    section.data(15).logicalSrcIdx = 14;
                    section.data(15).dtTransOffset = 14;

                    ;% rtB.iqcx4ji2gy
                    section.data(16).logicalSrcIdx = 15;
                    section.data(16).dtTransOffset = 15;

                    ;% rtB.b5jzfqhswx
                    section.data(17).logicalSrcIdx = 16;
                    section.data(17).dtTransOffset = 16;

                    ;% rtB.o5psjsgjbl
                    section.data(18).logicalSrcIdx = 17;
                    section.data(18).dtTransOffset = 17;

                    ;% rtB.plcvs3tnfq
                    section.data(19).logicalSrcIdx = 18;
                    section.data(19).dtTransOffset = 18;

                    ;% rtB.j0krl2hmjm
                    section.data(20).logicalSrcIdx = 19;
                    section.data(20).dtTransOffset = 19;

                    ;% rtB.l0b3kthfn5
                    section.data(21).logicalSrcIdx = 20;
                    section.data(21).dtTransOffset = 20;

                    ;% rtB.g5clq5nhx1
                    section.data(22).logicalSrcIdx = 21;
                    section.data(22).dtTransOffset = 21;

                    ;% rtB.mfk254enll
                    section.data(23).logicalSrcIdx = 22;
                    section.data(23).dtTransOffset = 22;

                    ;% rtB.id3dy3zvhz
                    section.data(24).logicalSrcIdx = 23;
                    section.data(24).dtTransOffset = 23;

                    ;% rtB.exiuxpqhzf
                    section.data(25).logicalSrcIdx = 24;
                    section.data(25).dtTransOffset = 24;

                    ;% rtB.jjqotddbo4
                    section.data(26).logicalSrcIdx = 25;
                    section.data(26).dtTransOffset = 25;

                    ;% rtB.hwx2zv05iw
                    section.data(27).logicalSrcIdx = 26;
                    section.data(27).dtTransOffset = 26;

                    ;% rtB.n3lttcqcta
                    section.data(28).logicalSrcIdx = 27;
                    section.data(28).dtTransOffset = 27;

                    ;% rtB.klsduu5tgr
                    section.data(29).logicalSrcIdx = 28;
                    section.data(29).dtTransOffset = 28;

                    ;% rtB.huebprkurf
                    section.data(30).logicalSrcIdx = 29;
                    section.data(30).dtTransOffset = 29;

                    ;% rtB.limtet2woe
                    section.data(31).logicalSrcIdx = 30;
                    section.data(31).dtTransOffset = 30;

                    ;% rtB.ap3x1johma
                    section.data(32).logicalSrcIdx = 31;
                    section.data(32).dtTransOffset = 31;

                    ;% rtB.iiwlxh5nit
                    section.data(33).logicalSrcIdx = 32;
                    section.data(33).dtTransOffset = 32;

                    ;% rtB.eiqmkhbtyt
                    section.data(34).logicalSrcIdx = 33;
                    section.data(34).dtTransOffset = 33;

                    ;% rtB.btfjf0eyh0
                    section.data(35).logicalSrcIdx = 34;
                    section.data(35).dtTransOffset = 34;

                    ;% rtB.cbizxpnxe4
                    section.data(36).logicalSrcIdx = 35;
                    section.data(36).dtTransOffset = 35;

                    ;% rtB.gzd0bvtkjr
                    section.data(37).logicalSrcIdx = 36;
                    section.data(37).dtTransOffset = 36;

                    ;% rtB.m3tgvpb0je
                    section.data(38).logicalSrcIdx = 37;
                    section.data(38).dtTransOffset = 37;

                    ;% rtB.kfmldcnjnq
                    section.data(39).logicalSrcIdx = 38;
                    section.data(39).dtTransOffset = 38;

                    ;% rtB.hdgz1v45gg
                    section.data(40).logicalSrcIdx = 39;
                    section.data(40).dtTransOffset = 39;

                    ;% rtB.lzxak0a220
                    section.data(41).logicalSrcIdx = 40;
                    section.data(41).dtTransOffset = 40;

                    ;% rtB.faofsw34yz
                    section.data(42).logicalSrcIdx = 41;
                    section.data(42).dtTransOffset = 41;

                    ;% rtB.oz3rry1khh
                    section.data(43).logicalSrcIdx = 42;
                    section.data(43).dtTransOffset = 42;

                    ;% rtB.jhxf0vx510
                    section.data(44).logicalSrcIdx = 43;
                    section.data(44).dtTransOffset = 43;

                    ;% rtB.ahjuj0uqhz
                    section.data(45).logicalSrcIdx = 44;
                    section.data(45).dtTransOffset = 44;

                    ;% rtB.eidzbo2axi
                    section.data(46).logicalSrcIdx = 45;
                    section.data(46).dtTransOffset = 45;

                    ;% rtB.dzr21vacs5
                    section.data(47).logicalSrcIdx = 46;
                    section.data(47).dtTransOffset = 46;

                    ;% rtB.g5dn334cgf
                    section.data(48).logicalSrcIdx = 47;
                    section.data(48).dtTransOffset = 47;

                    ;% rtB.hllpacqnya
                    section.data(49).logicalSrcIdx = 48;
                    section.data(49).dtTransOffset = 48;

                    ;% rtB.fad5yhgoij
                    section.data(50).logicalSrcIdx = 49;
                    section.data(50).dtTransOffset = 49;

                    ;% rtB.nzchyezcdl
                    section.data(51).logicalSrcIdx = 50;
                    section.data(51).dtTransOffset = 50;

                    ;% rtB.ouvddzwv14
                    section.data(52).logicalSrcIdx = 51;
                    section.data(52).dtTransOffset = 51;

                    ;% rtB.ixdckelvxz
                    section.data(53).logicalSrcIdx = 52;
                    section.data(53).dtTransOffset = 52;

                    ;% rtB.edzhheuutq
                    section.data(54).logicalSrcIdx = 53;
                    section.data(54).dtTransOffset = 53;

                    ;% rtB.ci3hl3yrta
                    section.data(55).logicalSrcIdx = 54;
                    section.data(55).dtTransOffset = 54;

                    ;% rtB.fgm5y4v5wu
                    section.data(56).logicalSrcIdx = 55;
                    section.data(56).dtTransOffset = 55;

                    ;% rtB.ciwvymzmbd
                    section.data(57).logicalSrcIdx = 56;
                    section.data(57).dtTransOffset = 56;

                    ;% rtB.bb2zjb32zc
                    section.data(58).logicalSrcIdx = 57;
                    section.data(58).dtTransOffset = 57;

                    ;% rtB.b0rncimetq
                    section.data(59).logicalSrcIdx = 58;
                    section.data(59).dtTransOffset = 58;

                    ;% rtB.eks3p0husl
                    section.data(60).logicalSrcIdx = 59;
                    section.data(60).dtTransOffset = 59;

                    ;% rtB.miudztndok
                    section.data(61).logicalSrcIdx = 60;
                    section.data(61).dtTransOffset = 60;

                    ;% rtB.aufytnaglh
                    section.data(62).logicalSrcIdx = 61;
                    section.data(62).dtTransOffset = 61;

                    ;% rtB.oguodtbfqn
                    section.data(63).logicalSrcIdx = 62;
                    section.data(63).dtTransOffset = 62;

                    ;% rtB.djfhnx1214
                    section.data(64).logicalSrcIdx = 63;
                    section.data(64).dtTransOffset = 63;

                    ;% rtB.afrqwefb1y
                    section.data(65).logicalSrcIdx = 64;
                    section.data(65).dtTransOffset = 64;

                    ;% rtB.n2srwnob2z
                    section.data(66).logicalSrcIdx = 65;
                    section.data(66).dtTransOffset = 65;

                    ;% rtB.folkmnz0t4
                    section.data(67).logicalSrcIdx = 66;
                    section.data(67).dtTransOffset = 66;

                    ;% rtB.k35lmpwnot
                    section.data(68).logicalSrcIdx = 67;
                    section.data(68).dtTransOffset = 67;

                    ;% rtB.gnjhdq1qvh
                    section.data(69).logicalSrcIdx = 68;
                    section.data(69).dtTransOffset = 68;

                    ;% rtB.hmcqeqif2o
                    section.data(70).logicalSrcIdx = 69;
                    section.data(70).dtTransOffset = 69;

                    ;% rtB.owtjgmweb5
                    section.data(71).logicalSrcIdx = 70;
                    section.data(71).dtTransOffset = 70;

                    ;% rtB.jpwjwhvvnr
                    section.data(72).logicalSrcIdx = 71;
                    section.data(72).dtTransOffset = 71;

                    ;% rtB.pju2ibjy53
                    section.data(73).logicalSrcIdx = 72;
                    section.data(73).dtTransOffset = 72;

                    ;% rtB.m0rkotwez0
                    section.data(74).logicalSrcIdx = 73;
                    section.data(74).dtTransOffset = 73;

                    ;% rtB.fsfzk25tf0
                    section.data(75).logicalSrcIdx = 74;
                    section.data(75).dtTransOffset = 74;

                    ;% rtB.jkrmpymttx
                    section.data(76).logicalSrcIdx = 75;
                    section.data(76).dtTransOffset = 75;

                    ;% rtB.mrxnnzezb3
                    section.data(77).logicalSrcIdx = 76;
                    section.data(77).dtTransOffset = 76;

                    ;% rtB.eqrybinqmz
                    section.data(78).logicalSrcIdx = 77;
                    section.data(78).dtTransOffset = 77;

                    ;% rtB.hx5qshmezi
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

                    ;% rtDW.j0cswhuhs0
                    section.data(1).logicalSrcIdx = 0;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            dworkMap.sections(1) = section;
            clear section

            section.nData     = 20;
            section.data(20)  = dumData; %prealloc

                    ;% rtDW.mg3u3mjjuf.LoggedData
                    section.data(1).logicalSrcIdx = 1;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.iljidhxeal.LoggedData
                    section.data(2).logicalSrcIdx = 2;
                    section.data(2).dtTransOffset = 2;

                    ;% rtDW.gfpi0a3qtl.LoggedData
                    section.data(3).logicalSrcIdx = 3;
                    section.data(3).dtTransOffset = 3;

                    ;% rtDW.p3eoc0sotd.LoggedData
                    section.data(4).logicalSrcIdx = 4;
                    section.data(4).dtTransOffset = 7;

                    ;% rtDW.naarlwo2lz.LoggedData
                    section.data(5).logicalSrcIdx = 5;
                    section.data(5).dtTransOffset = 8;

                    ;% rtDW.km0sbvoafb.LoggedData
                    section.data(6).logicalSrcIdx = 6;
                    section.data(6).dtTransOffset = 11;

                    ;% rtDW.lrvek4v1dp.LoggedData
                    section.data(7).logicalSrcIdx = 7;
                    section.data(7).dtTransOffset = 12;

                    ;% rtDW.homjuevopz.AQHandles
                    section.data(8).logicalSrcIdx = 8;
                    section.data(8).dtTransOffset = 13;

                    ;% rtDW.h1nrd1cc1r.AQHandles
                    section.data(9).logicalSrcIdx = 9;
                    section.data(9).dtTransOffset = 14;

                    ;% rtDW.l444m1qlra.AQHandles
                    section.data(10).logicalSrcIdx = 10;
                    section.data(10).dtTransOffset = 15;

                    ;% rtDW.fcjl52ojvi.AQHandles
                    section.data(11).logicalSrcIdx = 11;
                    section.data(11).dtTransOffset = 16;

                    ;% rtDW.cpiqtbmrnh.AQHandles
                    section.data(12).logicalSrcIdx = 12;
                    section.data(12).dtTransOffset = 17;

                    ;% rtDW.bya1bx33so.AQHandles
                    section.data(13).logicalSrcIdx = 13;
                    section.data(13).dtTransOffset = 18;

                    ;% rtDW.h1qcylur1s.AQHandles
                    section.data(14).logicalSrcIdx = 14;
                    section.data(14).dtTransOffset = 19;

                    ;% rtDW.ohqxmkyead.AQHandles
                    section.data(15).logicalSrcIdx = 15;
                    section.data(15).dtTransOffset = 20;

                    ;% rtDW.jnhubbfo5j.AQHandles
                    section.data(16).logicalSrcIdx = 16;
                    section.data(16).dtTransOffset = 21;

                    ;% rtDW.fokz15rol1.AQHandles
                    section.data(17).logicalSrcIdx = 17;
                    section.data(17).dtTransOffset = 22;

                    ;% rtDW.lj0dkmhahi.AQHandles
                    section.data(18).logicalSrcIdx = 18;
                    section.data(18).dtTransOffset = 23;

                    ;% rtDW.oqdb4h0bik.AQHandles
                    section.data(19).logicalSrcIdx = 19;
                    section.data(19).dtTransOffset = 24;

                    ;% rtDW.gkxaf1jyuy.LoggedData
                    section.data(20).logicalSrcIdx = 20;
                    section.data(20).dtTransOffset = 25;

            nTotData = nTotData + section.nData;
            dworkMap.sections(2) = section;
            clear section

            section.nData     = 16;
            section.data(16)  = dumData; %prealloc

                    ;% rtDW.atxpvhc5n0
                    section.data(1).logicalSrcIdx = 21;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.inlqdh4kh1
                    section.data(2).logicalSrcIdx = 22;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.glwrm0q1go
                    section.data(3).logicalSrcIdx = 23;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.apkag20puf
                    section.data(4).logicalSrcIdx = 24;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.mpaqmqlwie
                    section.data(5).logicalSrcIdx = 25;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.a4pqll03x1
                    section.data(6).logicalSrcIdx = 26;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.p1ctkshy11
                    section.data(7).logicalSrcIdx = 27;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.d34s0zbduf
                    section.data(8).logicalSrcIdx = 28;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.epcaebihmu
                    section.data(9).logicalSrcIdx = 29;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.kvbq5jwbvw
                    section.data(10).logicalSrcIdx = 30;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.audz3xgmq3
                    section.data(11).logicalSrcIdx = 31;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.loaiy2cacx
                    section.data(12).logicalSrcIdx = 32;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.e0qqpzfruh
                    section.data(13).logicalSrcIdx = 33;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.hamy5qs5qe
                    section.data(14).logicalSrcIdx = 34;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.je5l4cyavv
                    section.data(15).logicalSrcIdx = 35;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.ebwtwofxpd
                    section.data(16).logicalSrcIdx = 36;
                    section.data(16).dtTransOffset = 15;

            nTotData = nTotData + section.nData;
            dworkMap.sections(3) = section;
            clear section

            section.nData     = 10;
            section.data(10)  = dumData; %prealloc

                    ;% rtDW.aku4xs4fdc
                    section.data(1).logicalSrcIdx = 37;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.aigzt14viy
                    section.data(2).logicalSrcIdx = 38;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.jwuiy2atz3
                    section.data(3).logicalSrcIdx = 39;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.isc2qgvxqr
                    section.data(4).logicalSrcIdx = 40;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.awwuosv0rp
                    section.data(5).logicalSrcIdx = 41;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.g40nougyde
                    section.data(6).logicalSrcIdx = 42;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.cangmlit2r
                    section.data(7).logicalSrcIdx = 43;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.od3s0wtt1r
                    section.data(8).logicalSrcIdx = 44;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.ppnlaclhj1
                    section.data(9).logicalSrcIdx = 45;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.exujqdty5l
                    section.data(10).logicalSrcIdx = 46;
                    section.data(10).dtTransOffset = 9;

            nTotData = nTotData + section.nData;
            dworkMap.sections(4) = section;
            clear section

            section.nData     = 1;
            section.data(1)  = dumData; %prealloc

                    ;% rtDW.iln55qrgxk
                    section.data(1).logicalSrcIdx = 47;
                    section.data(1).dtTransOffset = 0;

            nTotData = nTotData + section.nData;
            dworkMap.sections(5) = section;
            clear section

            section.nData     = 16;
            section.data(16)  = dumData; %prealloc

                    ;% rtDW.i213auut4h
                    section.data(1).logicalSrcIdx = 48;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.e13oplclke
                    section.data(2).logicalSrcIdx = 49;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.laquyfl1kp
                    section.data(3).logicalSrcIdx = 50;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.f45mx2s03m
                    section.data(4).logicalSrcIdx = 51;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.gstdmwepbc
                    section.data(5).logicalSrcIdx = 52;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.ikvzyivi0j
                    section.data(6).logicalSrcIdx = 53;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.libi2dsvky
                    section.data(7).logicalSrcIdx = 54;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.ga2xty0o4o
                    section.data(8).logicalSrcIdx = 55;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.e1b1tx1d4l
                    section.data(9).logicalSrcIdx = 56;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.p0eq4zkt3w
                    section.data(10).logicalSrcIdx = 57;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.oahmxad2ft
                    section.data(11).logicalSrcIdx = 58;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.ksvy5jh2vc
                    section.data(12).logicalSrcIdx = 59;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.fyaulsukmr
                    section.data(13).logicalSrcIdx = 60;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.hdus1cnzep
                    section.data(14).logicalSrcIdx = 61;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.amiin4kdmj
                    section.data(15).logicalSrcIdx = 62;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.nwxy1ani3r
                    section.data(16).logicalSrcIdx = 63;
                    section.data(16).dtTransOffset = 15;

            nTotData = nTotData + section.nData;
            dworkMap.sections(6) = section;
            clear section

            section.nData     = 17;
            section.data(17)  = dumData; %prealloc

                    ;% rtDW.j5zzepflri
                    section.data(1).logicalSrcIdx = 64;
                    section.data(1).dtTransOffset = 0;

                    ;% rtDW.hk25kkxqwm
                    section.data(2).logicalSrcIdx = 65;
                    section.data(2).dtTransOffset = 1;

                    ;% rtDW.ffd0olavlx
                    section.data(3).logicalSrcIdx = 66;
                    section.data(3).dtTransOffset = 2;

                    ;% rtDW.eiohhcwuqa
                    section.data(4).logicalSrcIdx = 67;
                    section.data(4).dtTransOffset = 3;

                    ;% rtDW.eyyp3zs3le
                    section.data(5).logicalSrcIdx = 68;
                    section.data(5).dtTransOffset = 4;

                    ;% rtDW.nakeaw2rna
                    section.data(6).logicalSrcIdx = 69;
                    section.data(6).dtTransOffset = 5;

                    ;% rtDW.n4ne2vzdw3
                    section.data(7).logicalSrcIdx = 70;
                    section.data(7).dtTransOffset = 6;

                    ;% rtDW.ljfwvqlys1
                    section.data(8).logicalSrcIdx = 71;
                    section.data(8).dtTransOffset = 7;

                    ;% rtDW.ct4cduudvv
                    section.data(9).logicalSrcIdx = 72;
                    section.data(9).dtTransOffset = 8;

                    ;% rtDW.jkvzly4a1o
                    section.data(10).logicalSrcIdx = 73;
                    section.data(10).dtTransOffset = 9;

                    ;% rtDW.fxawe3htlm
                    section.data(11).logicalSrcIdx = 74;
                    section.data(11).dtTransOffset = 10;

                    ;% rtDW.irq32iym0t
                    section.data(12).logicalSrcIdx = 75;
                    section.data(12).dtTransOffset = 11;

                    ;% rtDW.jz3v2t25p2
                    section.data(13).logicalSrcIdx = 76;
                    section.data(13).dtTransOffset = 12;

                    ;% rtDW.m0i52xetkk
                    section.data(14).logicalSrcIdx = 77;
                    section.data(14).dtTransOffset = 13;

                    ;% rtDW.dfvcjzwy4h
                    section.data(15).logicalSrcIdx = 78;
                    section.data(15).dtTransOffset = 14;

                    ;% rtDW.nas3vtbw5t
                    section.data(16).logicalSrcIdx = 79;
                    section.data(16).dtTransOffset = 15;

                    ;% rtDW.awqg33it32
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


    targMap.checksum0 = 2453326653;
    targMap.checksum1 = 3278769589;
    targMap.checksum2 = 1795305681;
    targMap.checksum3 = 4235931226;

