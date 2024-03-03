ES -i genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 97 98 140
ES -i DadosZeEduardo_DATA.txt -m grn5ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 113 114 161
ES -i DadosZeEduardo_DATA.txt -m grn5ncyc -sd 0 113 114 161 -ind 0.509131859706313,4.705890160374611,0.100002011193567,0.100000000000000,0.587427439207285,0.930307173541365,0.271938661057075,0.995819658954900,0.440293368635028,0.294531357448824,0.643266490779396,0.100000206828619,0.999511289572582,0.443006784717132,0.999923201337715,0.999971822072094,0.230921764968310,0.999981525668626,1.027467752938993,20.590994565269891,1.000000000000000,24.932168591371969,1.000000000000000,1.809261972447725,1.603055004597870,25.000000000000000,24.960964294175983,12.005951470838969,25.000000000000000,15.878469782118251,24.523305696638271
ES -i GRN10_DATA.txt -m grn10new2 -e lsoda -a sade -n 1050000 -s 0 -sd 0 34 35 49

./batch_run.sh 4 30 "ES -i GRN10_DATA.txt -m grn10new2 -e lsoda -a sade -n 1050000 -s \$i -sd 0 34 35 49 >> grn10new2n0/output-grn10new2-sade.csv"
ES -i genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 97 98 140 -o output.txt -ind 3.550389036481703,4.999984009166924,2.362535860128819,0.209912729052144,0.660547030614924,0.332675428311095,0.477424728970194,0.999999999727546,0.202191953887012,0.571233265353244,0.579879784482317,0.885020452511245,0.540533078171088,0.559547118893801,0.100000235200446,2.078150429967071,21.763003547908450,22.834252918484779,24.937104685517361,6.653137007960251,17.257050315673268,2.949307820126937,24.760814468054135,24.041972510810019,18.478801009143232,1.124687462302994

ES -i GRN5_DATA.txt -m grn5new -e lsoda -a sade -n 1050000 -s 0 -sd 0 34 35 49 -o output.txt -ind 1.400038208318280,3.749405231944775,2.789412261097444,3.837757908866546,3.505168065738435,0.718175241938793,0.614596280026271,0.599593624645990,0.999998929122818,0.674924533676602,0.852759491371618,0.462658967629193,0.128494979617028,17.341925148523831,10.478756446847851,3.607855568756565,2.431778746682812,2.535361066863904,6.190938740899535,5.436562731634076,24.830196923712851

ES -i GRN10_DATA.txt -m grn10new2 -e lsoda -a sade -n 1050000 -s 0 -sd 0 34 0 34 -o output.txt -ind 0.612350570680964,1.773735401790358,1.491108196954083,2.898605139061402,0.140650910683491,3.530178368533732,0.328369686680488,0.101320991174217,2.008566425368409,0.416752065675240,0.674314752402210,0.567030086946192,0.935880424841606,0.423709533797724,0.963217790036262,0.788838719969977,0.112698347734706,0.513723653761847,0.981617594053271,0.531580757616678,0.483670957177960,0.733241669939404,0.248473162471431,0.154826933696811,0.133107494526149,0.357651106665407,0.522188692400949,0.987031032221768,0.496528298537652,0.387644920283008,0.750142770990936,0.948019045448747,0.107228182206860,0.205699739835332,0.682346493507898,1.827660700449575,12.817464860746856,4.502015725614627,7.971504198000627,22.929432872433245,1.157513860164513,23.817206903786480,24.799541799411805,19.620946189011196,7.221179829717094,2.271750143921959,18.209357921813410,9.137507002099265,24.660876963816772,14.862718789328825,23.413129412559808,3.837895703843942,1.047845896209167,2.459097084720614,3.328707050646093,16.042184760625563,24.064325745840466,11.078465701675607,2.047077655074208,16.751507287158322

ES -i genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 97 0 97 -o output.txt -ind 3.550389036481703,4.999984009166924,2.362535860128819,0.209912729052144,0.660547030614924,0.332675428311095,0.477424728970194,0.999999999727546,0.202191953887012,0.571233265353244,0.579879784482317,0.885020452511245,0.540533078171088,0.559547118893801,0.100000235200446,2.078150429967071,21.763003547908450,22.834252918484779,24.937104685517361,6.653137007960251,17.257050315673268,2.949307820126937,24.760814468054135,24.041972510810019,18.478801009143232,1.124687462302994

ES -i GRN5_DATA.txt -m grn5new -e lsoda -a sade -n 1050000 -s 0 -sd 0 34 0 34 -o output.txt -ind 1.400038208318280,3.749405231944775,2.789412261097444,3.837757908866546,3.505168065738435,0.718175241938793,0.614596280026271,0.599593624645990,0.999998929122818,0.674924533676602,0.852759491371618,0.462658967629193,0.128494979617028,17.341925148523831,10.478756446847851,3.607855568756565,2.431778746682812,2.535361066863904,6.190938740899535,5.436562731634076,24.830196923712851

ES -i DadosZeEduardo_DATA.txt -m grn5ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 113 0 113 -o output.txt -ind 3.244360864847833,3.898899828678533,0.101590149918321,0.100000000000719,1.012972162702005,0.977442198527650,0.307634599132096,0.454366619574674,0.155650064712432,0.781088230641081,0.649573154035218,0.999966819034412,0.605344220099708,0.999998028539363,0.999999999999559,0.892026773503344,0.341181356681746,0.112022829442903,1.045940688567378,1.821953501043235,24.250616671949853,21.442028538798549,16.714642667994372,18.648708445155143,24.007013077967105,1.209849245280861,6.778564462421043,24.681389365666568,9.735130528472181,1.797149257737955,24.480241194128112

ES -i DadosZeEduardo_DATA.txt -m grn5ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 113 0 113 -o output.txt -ind 0.259965643171702,0.541356773002716,2.723049408819733,0.100010503710478,0.683061143606153,0.100000000000000,0.896720597732577,0.100034088496540,0.998860923404447,0.100011911671525,0.644928372391290,0.103800078880908,0.999718441779722,0.987779922238183,0.999904434008379,0.956466062313451,0.131605754042801,0.212021585712222,25.000000000000000,1.101564794482040,24.999996457839909,24.983433967304279,4.719326554359199,1.043299836531002,13.831238950070755,1.317702340940321,24.999995778078805,25.000000000000000,3.421448278737768,1.564366329964603,2.894230361459349

ES -i genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 97 0 97 -o output.txt -ind 1.354396916894616,5.000000000000000,2.497461670696083,0.206547559117134,0.751815763421940,0.100000000000000,1.000000000000000,1.000000000000000,1.000000000000000,1.000000000000000,0.341683364407970,1.000000000000000,0.541740873257381,0.686739113180989,0.100000000000000,9.551987357138440,10.267541191171016,16.942559773670090,25.000000000000000,12.123108784867875,17.245641846654063,4.276603673496340,10.560481502493569,25.000000000000000,6.774313888927117,1.003622744752970

ES -i genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 97 98 140 -o output.txt -ind 1.354396916894616,5.000000000000000,2.497461670696083,0.206547559117134,0.751815763421940,0.100000000000000,1.000000000000000,1.000000000000000,1.000000000000000,1.000000000000000,0.341683364407970,1.000000000000000,0.541740873257381,0.686739113180989,0.100000000000000,9.551987357138440,10.267541191171016,16.942559773670090,25.000000000000000,12.123108784867875,17.245641846654063,4.276603673496340,10.560481502493569,25.000000000000000,6.774313888927117,1.003622744752970

ES -i genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 97 0 97 -o output.txt -ind 1.354396916894616,5.000000000000000,2.497461670696083,0.206547559117134,0.751815763421940,0.100000000000000,1.000000000000000,1.000000000000000,1.000000000000000,1.000000000000000,0.341683364407970,1.000000000000000,0.541740873257381,0.686739113180989,0.100000000000000,9.551987357138440,10.267541191171016,16.942559773670090,25.000000000000000,12.123108784867875,17.245641846654063,4.276603673496340,10.560481502493569,25.000000000000000,6.774313888927117,1.003622744752970

ES -i interpolated_train_genesABCD_data.txt -ts test_genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0

ES -i genesABCD_data.txt -m grn4ncyc -e lsoda -a cmaes -n 100000 -s 0 -sd 0 97 98 140 -o output.txt -ind 4.718765884081194,5.000000000000000,2.072607565237381,0.100000000000000,0.902633801007653,0.100000000000000,1.000000000000000,1.000000000000000,1.000000000000000,1.000000000000000,0.100000000000000,0.100000000000000,0.430199649078041,1.000000000000000,0.100000000000000,10.176538644901354,1.225758853022831,13.154692058139827,25.000000000000000,18.030823284144872,12.380684871718894,23.423754957181706,13.628870815782651,3.773602260062764,10.016183210347874,7.424283238894958

ES -i GRN5_DATA.txt -m grn5new -e lsoda -a sade -n 1050000 -s 0 -sd 0 34 35 49 -ro res_test.txt -po prog_test.txt