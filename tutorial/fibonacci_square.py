

from field import FieldElement
from polynomial import Polynomial, interpolate_poly, X, prod
from hashlib import sha256
from channel import serialize
from merkle import MerkleTree, verify_decommitment
from channel import Channel
from proof import Proof, DecommitmentData

def build_trace_eval():
  a = [FieldElement(1), FieldElement(3141592)]
  for _ in range(1021):
      v = a[-2]**2 + a[-1]**2
      a.append(v)
  assert a[-1] == FieldElement(2338775057)
  assert len(a) == 1023
  print("success generate trace")
  return a


# build trace domain
def build_domain(n):
  assert n % 2 == 0  # n must be even
  g = FieldElement.generator() ** ((FieldElement.k_modulus - 1) // n)
  assert g.is_order(n)

  # Fill G with the elements of G such that G[i] := g ** i
  G = [g ** i for i in range(n)]

  assert G[0] == FieldElement(1)
  assert G[1] == g
  for i in range(n//2):
    assert g * G[i] * G[-i-1] == FieldElement(1)

  print("success build trace domain")
  return G


def build_trace_domain():
  n = 1024
  G = build_domain(n)
  return G


def build_polynomial(xs, ys):
  f = interpolate_poly(xs, ys)
  v = f(2)
  assert v == FieldElement(1302089273)
  print("success build polynomial")
  return f


# build the larger evaluation domain x_coordinates
def build_lde_domain():
  n = 8192
  H = build_domain(8192)
  lde_domain = [H[i] * FieldElement.generator() for i in range(n)]

  assert len(set(lde_domain)) == len(lde_domain)
  w = FieldElement.generator()
  w_inv = w.inverse()
  assert '55fe9505f35b6d77660537f6541d441ec1bd919d03901210384c6aa1da2682ce' == sha256(
    str(H[1]).encode()).hexdigest(), \
    'H list is incorrect. H[1] should be h (i.e., the generator of H).'
  for i in range(n):
    assert ((w_inv * lde_domain[1]) ** i) * w == lde_domain[i]
  print('success build lde domain!')
  return lde_domain

def build_lde_eval(f, lde_domain):
  lde_eval = [f.eval(x) for x in lde_domain]

  assert '1d357f674c27194715d1440f6a166e30855550cb8cb8efeb72827f6a1bf9b5bb' == sha256(
    serialize(lde_eval).encode()).hexdigest()
  print('success build lde evaluation!')
  return lde_eval

def commit(data):
    tree = MerkleTree(data)
    return tree

def constraint_0(f):
    #constraint 1: f(x) = 1
    assert f.eval(1) == 1
    numerator = f - 1
    denominator = X - 1
    p0 = numerator / 	denominator
    r0 = numerator % 	denominator

    assert r0 == 0
    assert p0.eval(2718) == 2509888982
    assert p0.degree() == 1021
    print('Success build p0!')
    return p0

def constraint_1(f, x):
  assert f.eval(x) == 2338775057
  numerator = f - 2338775057
  denominator = X - x
  p1 = numerator / 	denominator
  r1 = numerator % 	denominator

  assert r1 == 0
  assert p1.eval(5772) == 232961446
  assert p1.degree() == 1021
  print('Success build p1!')
  return p1

def constraint_2(f,G):
  g = G[1]
  numerator = f.compose(g ** 2 * X) - f.compose(g * X) ** 2 - f ** 2
  denominator = (X**1024-1)/((X-g**1021)*(X-g**1022)*(X-g**1023))
  p2 = numerator / denominator
  r2 = numerator % denominator

  assert r2 ==0
  assert numerator(g ** 1020) == 0
  assert numerator.eval(g ** 1021) != 0

  assert p2.degree() == 1023, f'The degree of the third constraint is {p2.degree()} when it should be 1023.'
  assert p2.eval(31415) == 2090051528
  assert p2.degree() == 1023
  print('Success build p2!')
  return p2



def build_compostion_polynomial(p0, p1, p2, alpha0, alpha1, alpha2):
  # print('alpha0 =', alpha0, 'alpha1 =', alpha1, 'alpha2 =', alpha2)
  # print('build_compostion_polynomial, p0 =', p0.poly)
  # print('build_compostion_polynomial, p1 =', p1.poly)
  # print('build_compostion_polynomial,p2 =', p2.poly)
  cp = p0 * alpha0 + p1 * alpha1 + p2 * alpha2
  return cp

def next_fri_domain(domain):
  return [x**2 for x in domain[:len(domain)//2]]

def next_fri_polynomial(poly, beta):
  odd_coeffs = poly.poly[1::2]
  even_coeffs = poly.poly[::2]
  odd = beta * Polynomial(odd_coeffs)
  even = Polynomial(even_coeffs)
  return odd + even

def next_fri_layer(poly, domain, beta):
  next_domain = next_fri_domain(domain)
  next_poly = next_fri_polynomial(poly, beta)
  next_eval = [next_poly(x) for x in next_domain]
  return next_poly, next_domain, next_eval

def fri_commit(cp, lde_domain, cp_eval, cp_merkle, channel):
  fri_polys = [cp]
  fri_domains = [lde_domain]
  fri_layers = [cp_eval]
  fri_merkles = [cp_merkle]
  while fri_polys[-1].degree() > 0:
    beta = channel.receive_random_field_element()
    next_poly, next_domain, next_eval = next_fri_layer(fri_polys[-1], fri_domains[-1], beta)
    fri_polys.append(next_poly)
    fri_domains.append(next_domain)
    fri_layers.append(next_eval)
    fri_merkles.append(commit(next_eval))
    channel.send(fri_merkles[-1].root)  # record each layer's merkle root
  channel.send(str(fri_polys[-1].poly[0]))
  return fri_polys, fri_domains, fri_layers, fri_merkles

def decommit_on_fri_layers(fri_layers, fri_merkles, idx, channel):
  cp_proof = []
  for layer, merkle in zip(fri_layers[:-1], fri_merkles[:-1]):
    length = len(layer)
    idx = idx % length
    sib_idx = (idx + length // 2) % length
    channel.send(str(layer[idx]))
    channel.send(str(merkle.get_authentication_path(idx)))
    proof = DecommitmentData(idx, layer[idx], merkle.get_authentication_path(idx), merkle.root)
    cp_proof.append(proof)
    # valid = verify_decommitment(idx, layer[idx], merkle.get_authentication_path(idx), merkle.root)
    # assert valid == True, f'Failed to verify decommitment for f(g^i) at index {idx}.'
    channel.send(str(layer[sib_idx]))
    channel.send(str(merkle.get_authentication_path(sib_idx)))
    proof = DecommitmentData(sib_idx, layer[sib_idx], merkle.get_authentication_path(sib_idx), merkle.root)
    cp_proof.append(proof)
    # valid = verify_decommitment(sib_idx, layer[sib_idx], merkle.get_authentication_path(sib_idx), merkle.root)
    # assert valid == True, f'Failed to verify decommitment for f(g^i) at index {sib_idx}.'
  channel.send(str(fri_layers[-1][0]))

  return cp_proof, fri_layers[-1][0]


def decommit_on_query(lde_eval, tree, fri_layers, fri_merkles, idx, channel):
  assert idx + 16 < len(lde_eval), f'query index: {idx} is out of range. Length of layer: {len(lde_eval)}.'

  channel.send(str(lde_eval[idx]))  # f(x).
  channel.send(str(tree.get_authentication_path(idx)))  # auth path for f(x).
  x_proof = DecommitmentData(idx, lde_eval[idx], tree.get_authentication_path(idx), tree.root)
  # valid = verify_decommitment(idx, lde_eval[idx], tree.get_authentication_path(idx), tree.root)
  # assert valid == True, f'Failed to verify decommitment for f(x) at index {idx}.'

  next_idx = idx + 8
  channel.send(str(lde_eval[next_idx]))  # f(gx).
  channel.send(str(tree.get_authentication_path(next_idx)))  # auth path for f(gx).
  gx_proof = DecommitmentData(next_idx, lde_eval[next_idx], tree.get_authentication_path(next_idx), tree.root)
  # valid =verify_decommitment(next_idx, lde_eval[next_idx], tree.get_authentication_path(next_idx), tree.root)
  # assert valid == True, f'Failed to verify decommitment for f(gx) at index {next_idx}.'

  next_next_idx = idx + 16
  channel.send(str(lde_eval[next_next_idx]))  # f(g^2x).
  channel.send(str(tree.get_authentication_path(next_next_idx)))  # auth path for f(g^2x).
  g2x_proof = DecommitmentData(next_next_idx, lde_eval[next_next_idx], tree.get_authentication_path(next_next_idx), tree.root)
  # verify_decommitment(next_next_idx, lde_eval[next_next_idx], tree.get_authentication_path(next_next_idx), tree.root)
  # assert valid == True, f'Failed to verify decommitment for f(g^2x) at index {next_next_idx}.'

  cp_proof, final_eval = decommit_on_fri_layers(fri_layers, fri_merkles, idx, channel)

  proof = Proof(x_proof, gx_proof, g2x_proof, cp_proof, final_eval)
  return proof

def decommit_fri(channel):
  for query in range(3):
    # Get a random index from the verifier and send the corresponding decommitment.
    decommit_on_query(channel.receive_random_int(0, 8191 - 16), channel)

def test_next_fri_layer():
  test_poly = Polynomial([FieldElement(2), FieldElement(3), FieldElement(0), FieldElement(1)])
  test_domain = [FieldElement(3), FieldElement(5)]
  beta = FieldElement(7)
  next_p, next_d, next_l = next_fri_layer(test_poly, test_domain, beta)
  assert next_p.poly == [FieldElement(23), FieldElement(7)]
  assert next_d == [FieldElement(9)]
  assert next_l == [FieldElement(86)]
  print('success verify next_fri_layer!')



def original_stark101():
  print('original stark101 start!')
  trace_eval = build_trace_eval()  # [y0,y1,y2,...,y1022]
  trace_domain = build_trace_domain()  # [g^0,g^1,g^2,...,g^1023]

  f = build_polynomial(trace_domain[:-1], trace_eval)
  lde_domain = build_lde_domain()
  lde_eval = build_lde_eval(f, lde_domain)
  tree = commit(lde_eval)
  assert tree.root == '6c266a104eeaceae93c14ad799ce595ec8c2764359d7ad1b4b7c57a4da52be04'
  print('success commit lde evaluation!')

  channel = Channel()
  channel.send(tree.root)

  p0 = constraint_0(f)
  p1 = constraint_1(f, trace_domain[-2])
  p2 = constraint_2(f, trace_domain)

  #s0=b'-1021728838,1199291108,-459632242,-41103490,-431460358,1067145585,-1231134389,-72634232,1495414086,131992078,-1159696350,-939719258,182395720,-563918220,787363203,-1387478521,-317575747,940580168,-1452220647,172261911,194552229,-574211957,-1127044053,-1147687750,665276787,-1117791905,894601126,-1386121331,1414842856,1853390,1529169628,-365494408,-1363598387,-1500278976,1353072949,-4313449,-1193212177,-808896584,1449212945,-663885214,-1445110079,-240308321,-844489856,-130702230,-170244895,-114254846,-1422953326,-1051837630,1504084800,1403885957,1106418425,419968734,-1187420930,-496749109,-1330035703,-983073829,-1311647516,281065899,-697736463,-751430284,443188362,319061383,-1460026991,-1309151361,-1426088879,932125171,1392308025,1235404328,-344594224,47543930,-972176042,426657007,-1449883392,-783360991,1544322569,-777830020,-286018686,-370176170,-1279403703,-1451269284,1095213765,-41148398,189693363,-158034349,1070098519,1428544537,-1610577303,-500034224,-507383022,1063216470,-192283431,-29668822,30023706,471306309,-512119429,940164824,1026718907,384253082,-541590100,1472447937,559349125,-692101641,-161400927,1091303229,595768511,-1466924220,-1273089282,748330012,-1237919584,-1473628696,-589580291,-1448194544,1574501324,-1396776237,-210191671,1210328519,1171048400,784284975,583743363,-398322439,334047352,-88801874,911733644,317945957,1554844523,1493575160,1087647488,1048183580,-235588851,1043164579,225849423,-108999917,-1459560961,-192869852,1112463013,-1358564638,1054541359,-111936914,1342148998,1607997845,1452984441,-549416970,424384270,-1094973094,944391236,-680220850,-1489809829,1138224866,1490553160,675075551,325460109,686061647,751458096,1006715160,-1368082587,744562959,966988487,1218105406,-1610481173,988031589,-1182313277,1606574673,1373693099,-300328763,-105851492,837800003,-736684840,-512263722,-1215236476,1072199847,1069074664,884604808,-922769305,1042451657,89277372,356293629,-400998936,468100740,-473644982,1575220718,120272556,1590035980,653203856,1067502884,1005597702,534109458,1171648681,420653244,-334980661,-946459045,1331436439,-862196262,1571779964,-1080498922,1029773212,19028865,79827958,899026354,522991639,1399027273,301517489,1274578594,227118936,723518555,-137295309,-928880800,1399321197,-1190383445,330195153,807794577,1277389095,516910122,1059227677,-266404405,971879523,-1088555411,709593255,-1243083625,1592886473,-1195794402,-1608972733,602404994,-697018974,-365624957,-1496620821,-11251795,-64707293,920421805,1118342576,534749131,1143177326,-569006442,931251343,-14157164,585276201,1446578366,1253686226,1408975273,374666625,1129517923,941532004,328439626,526729975,1585790685,363607687,511731394,774624423,-1096842892,404106344,-133808164,1148156951,1550262977,-972906773,-96867432,-1346996508,1091273665,437076260,-540681857,-1223028619,-391315231,1501121886,898538461,-752681239,-292710225,1132270196,-1174694524,-891405801,-492615822,1442576735,1390640589,-954971912,-1274949772,450409893,-977009457,905501845,1264028093,569338553,-223568701,-1237382503,974191235,-625201026,-35336472,-440397964,731678076,-50058339,-157729648,451940606,-274782751,726824933,-98574944,916916156,1245881592,-190459809,-470527865,-164690884,1397091006,1307320515,-911096263,-537421787,-218486112,-1424423332,914280988,295523864,-59922270,178625772,123478856,-664475452,536519054,1429099107,-263749503,329669455,150776658,-1485755142,-287927255,758770131,-612782861,880812617,953785492,-105086853,1468831279,-492132518,-78738314,-820395206,-563953499,632254793,-530448049,-868485518,-283289757,709240827,-609809252,-306706702,-1498254222,-544122835,1241070474,-1138395005,-989898108,826884733,-984053338,-151419048,-595964290,-904930918,-1433972799,870525908,872337827,523160777,-606723566,-155243822,1213921383,-1360107612,-1181883571,-947003238,-729910952,-1313444539,-571905989,-577104497,198089653,555423142,-594164223,-963373703,1115130446,-612056951,798901174,669334281,-1486296921,1052245829,-1298014835,-739037686,-13402875,-256431103,-1160151325,1231871835,1028357488,-1967347,-1319416826,375031812,-1122156141,-1464696145,-946271445,1462191289,1455190670,-1061778539,373832353,-446962224,-181672260,760862139,-1495755608,872734679,164732371,133947145,-248286152,-1246286336,-1098810178,-987937523,677907634,-1547828438,108661396,-1185286963,-474231701,796291367,-765922579,971446437,-2620137,854467794,184852270,352264639,1231227708,602311157,1185853745,-1296371039,274845121,659789523,-795411126,450570192,949668391,602449779,1498611223,1184589682,-1127557360,583388998,-24033303,-1535291879,7386594,-120334790,-147562467,-467965278,143303981,501880707,-824386900,-880167207,906936898,646526432,1047271531,-1061377261,237411178,-1147293525,858771405,1582423545,-1524725516,-222439179,977462698,-1047129008,655415956,923148423,-1120429340,1137666985,323590362,1429570923,-903381308,502115844,1312659983,-430512070,-675576269,1058845562,-1440765551,755164113,623299152,-915605777,297937320,1280283168,229429089,-1206785441,1543899801,601867877,856845449,-1381893905,1125339057,-297735704,550782938,604047663,-870654836,-903590026,-1037096956,425926603,1580370515,177829935,-482972140,867869107,1406071342,-382148031,562572715,-450687914,-113922871,-606378505,-1531394262,1304402609,-1563466530,-558374953,-1576267608,-1565363148,-1533306645,1333473477,-904491909,-1293792894,-411036656,858403191,30855494,982510772,602174869,-1523247000,1609151159,-1301109444,435487918,400968233,-243425680,-1382364545,920130702,-617902269,-589885635,1491496635,94277547,1412552787,-461090091,1332887601,5948355,516516105,-304995888,-745648957,-1504930151,-589533127,-1395995084,-982429717,311822576,89093431,959739253,1370905163,-697056624,-788643805,1030497390,-110792236,-1241227158,933231155,675015489,517053615,972069872,1174054751,1109427112,-443033208,-1303769062,-3471366,1503298825,1433568208,28358320,-336207202,167848024,-1341176278,-298901391,-61955805,208828338,1024372629,-407998323,306855290,-1244814107,-259003461,-535676781,867234552,-21320428,-491354903,-1565192718,-258659667,-297599942,336233860,222290189,-917074531,1451552785,826170630,-1168916725,-1324759243,421588123,-1166164487,-1310229661,-145236940,1380924690,902197161,752412379,436111697,-509883809,21188312,243094457,1606618970,-878815508,249985814,-256593391,-621133769,1595817631,167980028,1134437377,1158493557,-856104932,1303179627,121061321,380820127,-1547356217,879102415,960025543,953984174,-544049567,471341325,-38700934,822674737,247125194,-558636541,-8948925,938952390,-925006622,1104259059,1144130116,1034371860,1173738351,1408651578,563935237,1181529411,1002942135,1281745288,957208115,1013069972,50709163,-1146700663,1156441238,62647318,-1018843361,-233322439,1250205151,946183037,249516951,-1151795731,-428098209,-1354132779,-825754300,11533624,-825771347,1287385692,-800557214,1509147792,-1144068488,556258318,1489695223,-180042777,-833476737,-61370138,-1433458865,1278485157,1007474982,1570247456,-1285433757,-198728459,357576888,1456328852,-324087677,-1524660931,-260170169,-1288132829,802442645,-938409717,35407348,1463878911,-785409189,-31192848,-1453439030,-620199683,-872918459,-968036297,-863214396,1571564564,1465680812,483255402,909952486,-70882741,-1562567887,-772134742,-779139042,-100828904,726454124,811999293,747240062,1590494480,736714099,1208438954,542827979,-691615565,1296952745,-708755444,-811235899,49551866,1084628928,-1140443246,-596108598,690772622,430207157,807021485,425240026,-107345204,126516163,-1419604091,-179997574,-342921885,149503113,1076808629,170337536,-851267924,-665315581,837125888,138275212,-1089964298,87456016,-823792456,1339418025,-1211452076,-476934506,-734051746,-635175722,-687569429,609711844,-1079949200,1033418148,675875289,1038063283,-105658954,-234766270,-772625280,-783428922,798158970,500487024,467423079,1041548272,-754502566,-690973528,-147811552,609262273,1426684548,902264566,1344827832,1212945371,1007033793,519975270,-1187000185,611328835,1206699673,-236384422,54949139,126887530,-1084749431,960623871,-1389192516,862286121,1021013986,746885664,-971659115,733006448,963593796,905666571,259907397,-1398807532,-45168674,636240857,-1016142799,-789988725,-610235076,35260429,369460637,797400544,-1094680246,-2344467,-268447079,-1034412356,-820966392,465435949,1608019351,-763918560,-1176457242,807800538,-1260311657,1138686919,638496662,224777363,-1121361119,1048950245,-463718320,-62115865,1017374826,-718909122,-447353684,-1006989138,-957670233,54963754,-785755649,-1153839799,202387494,-1236365200,571833505,131537225,-943250897,-253718640,1006450175,-67388722,4994287,939978680,1480332076,1544372737,-1346356515,-793834081,614827184,229976707,859683791,-60473257,1312774031,-1138347450,1025542896,-1146882936,-700480936,-95028057,1173075858,328744862,-218359253,-600942341,325925743,-447348269,-1198995307,385852610,-92644251,-70896099,943714212,753565439,57282309,211198498,40326536,1367382363,1100340381,-7952874,-140213158,112259150,-615397775,504928132,948923260,1571798608,748078185,1205160946,-1161541491,-1420997801,57701363,1010992144,1084772250,-1368471009,-468972138,-1522420239,-335574118,853071038,-448007570,989286771,-1264271876,1529010769,-894001733,-222070000,1201975559,920100518,1445926135,-582225590,-655104066,-1594266731,541708257,-1506278007,-998492987,-1052843757,-392967103,1507494586,1099725981,90560388,638975241,950946989,248751631,-369940450,-787652901,-1328459495,-1514826263,1144234014,-284776854,1270730227,-1085788077,-363604168,-770348390,316530348,-511289166,-1191677878,934587713,1430850762,843353984,234827380,1519220898,-423399563,558617,1294152769,-790246814,1096882338,811673643,263321422,-1498338557,365966302,-354290053,-1204345382,367075086,-48943879,-106439658,-1329305985,-1421099786,-136312293,197345272,955985287,360244106,-546566314,-395440856,756864257,-279744929,286070524,-1055505831,1099925804,185136828,-38749668,866902134,417434769,-25805891,1142745924,212628441,435762059,-522163551,-639393634,1234182006,-747518115,1145817385,918584374,30779341,1263531399,1428410855,1028370183,-641850399,1387837289,-316556697,737337679,-1407471277,-1429917823,1331953756,-345128126,-981890962,709322567,799351072,-113492937,827973287,1233185200,458621149,-875788332,-1530341014,1299898849,418224582,-1475976659,1187486926,619978774,1553921853,-890074822,-508882038,-993648893,877751494,-1526319694,-126592598,648604100,1085342665,1126957631,987718322,1213121929,-815870268,-1171647912,-331895053,843863494,941298087,1264845884,1145992411,462905904,1059650771,1335575459,911861158,88378396,-793376551,534451302,754905623,4144858,-1606601804,-315388973,213857542,1406958956,-1120990766,-1303359428,1176255667,735073659,1592559896,-425977635,1383991623,923246118,-1470707321,521708992,-325116392,-836074820,-1550598713,770704514,-974827393,29420497,218060227,-459546018,-1440716102,-603788727,306075179,-467731229,1249900805,484544273,-773626289,-699871262,778626990,1022485444,-617864856,364011434,1435216806,-215795471,615669792,1005459075,1127986164,934287394,-72906227,1145297574'
  alphas = channel.derive_alphas(3)
  cp = build_compostion_polynomial(p0, p1, p2, alphas[0], alphas[1], alphas[2])
  print('cp =', cp.poly)
  cp_eval = [cp(x) for x in lde_domain]

  cp_tree = commit(cp_eval)
  assert cp_tree.root == 'd7e5200e990727c6da6bf711aeb496244b8b48436bd6f29066e1ddb64e22605b'
  print('success commit cp evaluation!')

  channel.send(cp_tree.root)
  fri_polys, fri_domains, fri_layers, fri_merkles = fri_commit(cp, lde_domain, cp_eval, cp_tree, channel)
  # proof= ['send:6c266a104eeaceae93c14ad799ce595ec8c2764359d7ad1b4b7c57a4da52be04', 'send:d7e5200e990727c6da6bf711aeb496244b8b48436bd6f29066e1ddb64e22605b', 'send:dcbb68ffbf425707f19f44d6c7e48026bb08fc977cfc91ca82cc29b89913c037', 'send:b59994a75458d007e5b67f02700b7f96bbc510900c1fb6124797575318a7340c', 'send:92b9297e4d1920f6c91500350a22592e41f35d3e127bac0d68c36224cb3e8abe', 'send:43b2754838e9f452f7dd8fe0b7ebd642e1e157d27ca810775b19a49e74512434', 'send:c26be2d587fad765c7c81d888d5f9d8a3c92f8ae4ef802c4b1be6c6b483dba78', 'send:4c7e733cf2c6ed1541a0d3dd11f8df5a736d4fa5410114f9f7187253f0beaf3c', 'send:f42d6ba6d5a1e8ab9b8c2e44eb9fd8dd296c949a1fb8152fda8c9398bea2474f', 'send:10f032ad5a87c7a432b035b3b759bdef277a851b3ed99296cb89c73205a5a605', 'send:cd47fddeda8d081481286de849c04a4746f1541d3f047a7d56837ed5e517d7d8', 'send:57fbccfbb49bc6389f4327bcfc9ae7b32a9bd6593ddcc30238f80a399a45906c', 'send:807809296']

  proof_1 = decommit_on_query(lde_eval, tree, fri_layers, fri_merkles, 1, channel)
  proof_100 = decommit_on_query(lde_eval, tree, fri_layers, fri_merkles,100, channel)

  proof_1.verify(proof_1.final_eval)
  proof_100.verify(proof_1.final_eval)
  print('success original stark101!')


if __name__ == "__main__":
    original_stark101()



