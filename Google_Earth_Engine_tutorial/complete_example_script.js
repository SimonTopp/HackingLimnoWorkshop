// Google Earth Engine script prepared for AEMON-J Remote Sensing day workshop
// by Xiao Yang (yangxiao@live.unc.edu)
// license: MIT

// ******************************************************************************** //
// useful resources, 1hr is not enough to learn GEE but we will get a taste
// ******************************************************************************** //

// ******************************************************************************** //
// tour of the interface
// ******************************************************************************** //

// "Feature tour" function from Google Earth Engine

// ******************************************************************************** //
// WORKING WITH OPTICAL REMOTE SENISNG DATA
// ******************************************************************************** //
// data import (using the search bar, pick the collection, and import)
var landsat8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR");
var point = ee.Geometry.Point([-112.6674, 41.2205]);

print(landsat8.first()); // ee.Image consists of bands and properties

landsat8 = landsat8
.filterMetadata('CLOUD_COVER', 'less_than', 30)
.filterBounds(point);

var img = ee.Image(landsat8.first());

Map.addLayer(img, {bands: ['B4', 'B3', 'B2'], min: 0, max: 10000, gamma: 1.2}, 'rgb example', false);

// cloud and cloud shadow (make use the example cloud-masking scripts in the "Examples" repo)
// snow/ice (masking snow/ice in addition to clouds)

// Function to cloud mask from the pixel_qa band of Landsat 8 SR data.
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = 1 << 3;
  var snowIceBitMask = 1 << 4;
  var cloudsBitMask = 1 << 5;

  // Get the pixel QA band.
  var qa = image.select('pixel_qa');

  // get bands of cloud/shadow and ice for flagging the results
  var cloud = qa.bitwiseAnd(cloudShadowBitMask).neq(0).or(qa.bitwiseAnd(cloudsBitMask).neq(0)).rename('cloud');
  var snowIce = qa.bitwiseAnd(snowIceBitMask).neq(0).rename('snowIce');

  // Both cloud and snowIce flags should be set to zero, indicating clear conditions.
  var mask = cloud.not()
      .and(snowIce.not());

  // Return the masked image, scaled to reflectance, without the QA bands.
  return image.updateMask(mask).divide(10000)
      .select("B[0-9]*")
      .addBands(cloud)
      .addBands(snowIce)
      .copyProperties(image)
      .copyProperties(image, ['system:time_start']);
}

// water classification example
var Mndwi = function(image) {
  return(image.normalizedDifference(['B3', 'B6']).rename('mndwi'));
};

var imgMasked = ee.Image(maskL8sr(img));
print(imgMasked, 'imgMasked');
var waterMask = Mndwi(img).gte(0);
print(waterMask, 'waterMask');

Map.addLayer(waterMask.selfMask(), {palette: ['blue']}, 'water', false);
Map.addLayer(imgMasked.select(['cloud']).selfMask(), {palette: ['grey']}, 'cloud', false);
Map.addLayer(imgMasked.select(['snowIce']).selfMask(), {palette: ['cyan']}, 'snowIce', false);

// terrain shadow (particularly important for alpine lakes)

var CalcHillShadowSR = function(image) {
  var dem = ee.Image("MERIT/DEM/v1_0_3").select('dem');//.clip(image.geometry().buffer(9000).bounds());
  var SOLAR_AZIMUTH_ANGLE = ee.Number(image.get('SOLAR_AZIMUTH_ANGLE'));
  var SOLAR_ZENITH_ANGLE = ee.Number(image.get('SOLAR_ZENITH_ANGLE'));
  var hillShadow = ee.Terrain.hillShadow(dem, SOLAR_AZIMUTH_ANGLE, SOLAR_ZENITH_ANGLE, 100, true).rename(['hillshadow']);
  return(image.updateMask(hillShadow).addBands(hillShadow));
};

imgMasked = CalcHillShadowSR(imgMasked);// output 1 means surface illuminated.
Map.addLayer(imgMasked, {bands:['hillshadow'], min: 0, max: 1}, 'terrainShadow', false);


// ******************************************************************************** //
// EXAMPLE SCENARIO 1: reflectance time series over a sample location (color)
// ******************************************************************************** //

var rgb2wavelengthLehmann = function(img, rgbNames) {

  rgbNames = ee.List(rgbNames);

  var bandMap = {
    'U': img.select(rgbNames.getString(0)),
    'R': img.select(rgbNames.getString(1)),
    'G': img.select(rgbNames.getString(2)),
    'B': img.select(rgbNames.getString(3))
  };

  var xi = img.expression('34.457*R + 51.135*G + 6.950*B + 11.053*U', bandMap);
  var yi = img.expression('18.034*R + 66.023*G + 21.053*B + 1.320*U', bandMap);
  var zi = img.expression('0.016*R + 2.606*G + 34.931*B + 58.038*U', bandMap);
  var xyzsum = xi.add(yi).add(zi);
  var x = xi.divide(xyzsum).subtract(0.3333);
  var y = yi.divide(xyzsum).subtract(0.3333);
  var z = zi.divide(xyzsum);

  // test .atan function: https://code.earthengine.google.com/5de0fdaec32709ce4c48eb39505d3f7a
  var alpha = y.atan2(x).multiply(180).divide(Math.PI).rename('alpha'); //note the difference in input order between atan2 in R and in GEE

  // reverse the alpha so that it increase counterclockwise and start from x-positive
  var alphaRotated = alpha.multiply(-1).add(90);

  var alpha100 = alphaRotated.divide(100);
  var deltaAlpha = alpha100.expression('-52.16*(a**5) + 373.81*(a**4) - 981.83*(a**3) + 1134.19*(a**2) - 533.61*a + 76.72', {'a': alpha100.select('alpha')});
  var alphaRotatedCorrected = alphaRotated.add(deltaAlpha);

  // converted it back to increase clockwise and starting from y-negative
  alpha = alphaRotatedCorrected.multiply(-1).add(90);

  var xRef = [-154.377518550783, -154.374079185343, -154.371310981132, -154.368399454696, -154.365115251439, -154.361447077962, -154.357458377889, -154.35317248265, -154.348182658999, -154.34256905297, -154.33633712409, -154.328919388877, -154.321078124016, -154.313438423011, -154.306700024077, -154.301287233547, -154.296633126868, -154.292117937264, -154.287821123092, -154.28311210129, -154.278064907065, -154.27224522106, -154.265150147062, -154.256277138, -154.245912632646, -154.234411874566, -154.221568326131, -154.207879696225, -154.194188685055, -154.181273844653, -154.169347120844, -154.157064018022, -154.143724679095, -154.128974738246, -154.113311381481, -154.097093006604, -154.080600188505, -154.062422849999, -154.040430178147, -154.013562000451, -153.981531776792, -153.945327702279, -153.90643779259, -153.865356052419, -153.822081845801, -153.77618649188, -153.727522766462, -153.675730674166, -153.62023573598, -153.560322319785, -153.4956992656, -153.425436286933, -153.34980819458, -153.269514264977, -153.184896650184, -153.096151212412, -153.003118934376, -152.904997596748, -152.801421737405, -152.692378801833, -152.577292490058, -152.455940553747, -152.327799442904, -152.19155812278, -152.046546005902, -151.892078476801, -151.72733184776, -151.552192830568, -151.366698191189, -151.170808147938, -150.964777795169, -150.748926947012, -150.523180533515, -150.286294415701, -150.037203969858, -149.774181045867, -149.495584236907, -149.200866799776, -148.89023396176, -148.563725068977, -148.220956720414, -147.861488104698, -147.48217538719, -147.077402171674, -146.64044046865, -146.162116108496, -145.63163176621, -145.042885510952, -144.39190936266, -143.676567427747, -142.897184683129, -142.052478360742, -141.13251317322, -140.127175294255, -139.024552164264, -137.809739298933, -136.468711202927, -134.988903009653, -133.356499285147, -131.556577701595, -129.574799824384, -127.393772236356, -124.993574978312, -122.356122923578, -119.468124981564, -116.324219852632, -112.930378901234, -109.291290774868, -105.412505201446, -101.310758944573, -97.0155498546549, -92.5771777449705, -88.0701019901376, -83.5725721528921, -79.1607324654358, -74.9013659090351, -70.8494434767422, -67.0334847070785, -63.4641089663727, -60.1445504878684, -57.0723296661657, -54.2469054053975, -51.6622774128365, -49.2972443895531, -47.1249887368399, -45.1160060067983, -43.2469814955475, -41.5064414303078, -39.8839012823132, -38.3699626736652, -36.9562400333347, -35.6318522542182, -34.390274812277, -33.2290396542123, -32.144672368948, -31.132876380846, -30.1859738587168, -29.2955930809701, -28.4563721735574, -27.6628096672458, -26.9094192905555, -26.1912119514419, -25.5001355090307, -24.8274325436337, -24.1655666130937, -23.5079678040273, -22.8505302532628, -22.1925838152901, -21.5348240760989, -20.8773037299329, -20.2205852765341, -19.5638112768114, -18.9041602407473, -18.2397532144779, -17.5688036182667, -16.8893860178191, -16.2000120351211, -15.4993514035446, -14.7853476872828, -14.0560552966185, -13.3093189410279, -12.5431137313446, -11.7557451688982, -10.9457859001674, -10.1112084195601, -9.25074433495616, -8.3621370562045, -7.44333003671669, -6.49188345150389, -5.5053183851507, -4.48141320145393, -3.41714395762067, -2.31084559049961, -1.16025968431981, 0.0368738206174777, 1.28265444212602, 2.58043202010527, 3.93215895375249, 5.3396293956838, 6.80407318410594, 8.3278620377684, 9.9129316864974, 11.5592308319806, 13.2677745205127, 15.0386715586981, 16.8719840494074, 18.7670295925627, 20.7213840750234, 22.7326769510584, 24.7972811129161, 26.9117764484611, 29.0714674983749, 31.2702482728973, 33.5014999172376, 35.7581535377841, 38.0323298424252, 40.3157735164535, 42.5998562549555, 44.8756197609734, 47.133984147698, 49.3671492969771, 51.5680789211409, 53.729593604541, 55.8452011206872, 57.908598979884, 59.9147112053308, 61.8585768867874, 63.7378681602241, 65.5506232837678, 67.2962033411433, 68.9740636483888, 70.5852753052773, 72.1278393940881, 73.5973889029212, 74.9913612828952, 76.3073107167814, 77.5451600817744, 78.7113679253691, 79.8141974701489, 80.8612945218091, 81.8590960379573, 82.812150382453, 83.7197405183344, 84.5799285763649, 85.3918135292913, 86.1537575409584, 86.8667228656646, 87.5352988484769, 88.1637370444293, 88.7559910480277, 89.3155834432002, 89.8446186715729, 90.3440799314638, 90.8152031631223, 91.259298776048, 91.6775806245166, 92.0711784238747, 92.4417800069171, 92.7906496876107, 93.1191019955034, 93.4285388569152, 93.7201162619478, 93.9953985483001, 94.2557228762623, 94.502080117381, 94.7352623134217, 94.9560495107269, 95.1654532610661, 95.3646766021536, 95.5550372807154, 95.7377914040296, 95.9136461598127, 96.0826222684015, 96.2444265566319, 96.3990795449417, 96.5465490849891, 96.6874012171133, 96.8218445716112, 96.9499663273297, 97.0720286189983, 97.1882145488163, 97.2990054551434, 97.4044767071213, 97.5044310046874, 97.5988524945024, 97.6873997654724, 97.7696470213888, 97.8461739914784, 97.9175294273907, 97.984475961602, 98.0480873591, 98.1092635533621, 98.1674805814859, 98.2227616488076, 98.2748877023783, 98.3240014108078, 98.3701230276047, 98.4137504655831, 98.4545401856561, 98.4923889268348, 98.5270739830665, 98.5586106870575, 98.5871317149722, 98.6131236062136, 98.6371880008124, 98.6600415258744, 98.6826338546069, 98.7047295874307, 98.7263302427477, 98.7470837283615, 98.7667570180146, 98.7852354678059, 98.8028748467829, 98.8196772474933, 98.8358793855893, 98.8520689968486, 98.8684804519402, 98.884996125723, 98.9016157533379, 98.9178714654897, 98.9337640936549, 98.9494111758912, 98.964813248644, 98.9799708395498, 98.9943021043201, 99.0072267336847, 99.018398321448, 99.0274708205407, 99.0350282137731, 99.0414207858229, 99.0471143305894, 99.0526901682721, 99.0578000384308, 99.0625603771204, 99.0668553559912, 99.0705692159784, 99.0739343278111, 99.0771828837435, 99.0799669556229, 99.0819387805879];
  var yRef = [380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600, 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 698, 699];

  // add purity
  var s = x.multiply(x).add(y.multiply(y)).sqrt();
  var deltaS = alpha100.expression('-0.0099*(a**5) + 0.1199*(a**4) - 0.4594*(a**3) + 0.7515*(a**2) - 0.5095*a + 0.1222', {'a': alpha100.select('alpha')});
  s = s.add(deltaS).rename('dist2wp');

  return(alpha.interpolate(xRef, yRef, "mask").rename('dwLehmann').addBands(s));
};
var colorVis = {
      min: 471,
      max: 600,
      palette: ['#2158bc', '#2158bc', '#2158bc', '#2158bc', '#2158bc', '#316dc5', '#316dc5', '#316dc5', '#316dc5', '#316dc5', '#327cbb', '#327cbb', '#327cbb', '#327cbb', '#327cbb', '#4b80a0', '#4b80a0', '#4b80a0', '#4b80a0', '#568f96', '#568f96', '#568f96', '#568f96', '#568f96', '#568f96', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#6d9298', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#698c86', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#759e72', '#7ba654', '#7ba654', '#7ba654', '#7ba654', '#7ba654', '#7ba654', '#7ba654', '#7ba654', '#7ba654', '#7ba654', '#7dae38', '#7dae38', '#7dae38', '#7dae38', '#7dae38', '#94b660', '#94b660', '#94b660', '#94b660', '#a5bc76', '#aab86d', '#adb55f', '#a8a965', '#a8a965', '#ae9f5c', '#ae9f5c', '#b3a053', '#b3a053', '#af8a44', '#af8a44', '#a46905', '#a46905', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04', '#9f4d04']
  };

var landsat8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR");

print(landsat8.first()); // ee.Image consists of bands and properties

landsat8 = landsat8
.filterMetadata('CLOUD_COVER', 'less_than', 30)
.filterBounds(poi)
.map(maskL8sr)
.map(CalcHillShadowSR)
.map(function(i) {
  return(i.addBands(rgb2wavelengthLehmann(i, ['B1', 'B4', 'B3', 'B2'])));
});

// on one image
var img = ee.Image(landsat8.first());
var waterMask = Mndwi(img).gte(0);

print(img, 'img', ee.Date(img.get('system:time_start')));
Map.addLayer(img, {bands: ['B4', 'B3', 'B2'], min: 0, max: 1, gamma: 1.8}, 'rgb example', false);
Map.addLayer(img.select(['dwLehmann']).updateMask(waterMask), colorVis, 'water color', false);

// get reflectances value via inspector

// extract reflectance value using reduceRegion

print(img.reduceRegion(ee.Reducer.first(), poi, 30));

// get time series using mapped function

var ExtractReflectance = function(image) {
  var dict = image.reduceRegion(ee.Reducer.first(), poi, 30);
  var doy = ee.Date(image.get('system:time_start')).format('D');
  return(ee.Feature(poi, dict).copyProperties(image, ['system:time_start']).set('doy', ee.Number.parse(doy)));
};

print(ExtractReflectance(img));

var targetBand = 'dwLehmann';
var refTs = landsat8.map(ExtractReflectance)
.filterMetadata(targetBand, 'not_less_than', 0);

var tsChart = ui.Chart.feature.byFeature(refTs.sort('system:time_start'), 'system:time_start', [targetBand]);
print(tsChart);
var doyChart = ui.Chart.feature.byFeature(refTs.sort('doy'), 'doy', [targetBand]).setChartType('ScatterChart');
print(doyChart);
