var plotsBaseURL = '../3-Var_modes/'
var imgbasename = '3m'

var origin = document.querySelector("#sel_origin");
var stdate = document.querySelector("#sel_stdate");
var varn = document.querySelector("#sel_varn");
var score = document.querySelector("#sel_score");
var season = document.querySelector("#sel_season");

var origin_valor = "ecmwf.s51"
var stdate_valor = "11"
var varn_valor = "nao_box"
var score_valor = "corr"
var folder2 = ""
var freq = "monthly"
var season_valor = "DJF"

const seas_array = ["JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ", "DJF"];
const month_array = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];


origin.addEventListener('input',function(evento){
    origin_valor = evento.target.value
    console.log(origin_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
        var folder = "1-Box_calc/"
      } else {
        var folder = "2b-EOFs_calc/"
    }

    // if (str_stdate === "08") {
    //     var folder2 = "daily_data_extension/"
    //     var freq = "daily"
    //   } else {
    //     var folder2 = ""
    //     var freq = "monthly"
    // }

    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.png'
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
    })
stdate.addEventListener('input',function(evento){
    stdate_valor = evento.target.value
    console.log(stdate_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
        var folder = "1-Box_calc/"
      } else {
        var folder = "2b-EOFs_calc/"
    }

    // if (str_stdate === "08") {
    //     var folder2 = "daily_data_extension/"
    //     var freq = "daily"
    //   } else {
    //     var folder2 = ""
    //     var freq = "monthly"
    // }
      
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.png'
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
  })
varn.addEventListener('input',function(evento){
    varn_valor = evento.target.value
    console.log(varn_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
        var folder = "1-Box_calc/"
        var eof_num = "1"
      } else if (varn_valor === "nao_eof") {
        var folder = "2b-EOFs_calc/"
        var eof_num = "1"
      } else if (varn_valor === "ea") {
        var folder = "2b-EOFs_calc/"
        var eof_num = "2"
      } else if (varn_valor === "eawr") {
        var folder = "2b-EOFs_calc/"
        var eof_num = "3"
      } else if (varn_valor === "sca") {
        var folder = "2b-EOFs_calc/"
        var eof_num = "4"
    }

    // if (str_stdate === "08") {
    //     var folder2 = "daily_data_extension/"
    //     var freq = "daily"
    //   } else {
    //     var folder2 = ""
    //     var freq = "monthly"
    // }

    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.png'
    const plot3 = document.querySelector("#plot3") 
    plot3.src = plotsBaseURL+'2-EOFs_calc/plots/ERA5_EOF'+eof_num+'_'+season_valor+'.png'  
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
})
score.addEventListener('input',function(evento){
    score_valor = evento.target.value
    console.log(score_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
        var folder = "1-Box_calc/"
      } else {
        var folder = "2b-EOFs_calc/"
    }

    // if (str_stdate === "08") {
    //     var folder2 = "daily_data_extension/"
    //     var freq = "daily"
    //   } else {
    //     var folder2 = ""
    //     var freq = "monthly"
    // }

    const plot1 = document.querySelector("#plot1") 
    plot1.src = plotsBaseURL+'3-Common_plots/Score-card_'+score_valor+'_'+season_valor+'_3m.png'
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
})

season.addEventListener('input',function(evento){
  season_valor = evento.target.value
  console.log(season_valor)
  var str_stdate = String(stdate_valor).padStart(2, '0');

  if (varn_valor === "nao_box") {
    var folder = "1-Box_calc/"
    } else {
    var folder = "2b-EOFs_calc/"
  }

  if (varn_valor === "nao_box") {
    var folder = "1-Box_calc/"
    var eof_num = "1"
  } else if (varn_valor === "nao_eof") {
    var folder = "2b-EOFs_calc/"
    var eof_num = "1"
  } else if (varn_valor === "ea") {
    var folder = "2b-EOFs_calc/"
    var eof_num = "2"
  } else if (varn_valor === "eawr") {
    var folder = "2b-EOFs_calc/"
    var eof_num = "3"
  } else if (varn_valor === "sca") {
    var folder = "2b-EOFs_calc/"
    var eof_num = "4"
  }

  const isthisSeason = (element) => element == season_valor;
  var season_index = seas_array.findIndex(isthisSeason)
  var lead1_index = season_index-1
  if(lead1_index<0){lead1_index = lead1_index+12}
  var lead2_index = season_index-2
  if(lead2_index<0){lead2_index = lead2_index+12}
  var lead3_index = season_index-3
  if(lead3_index<0){lead3_index = lead3_index+12}

  var lead1 = document.querySelector("#lead1")
  lead1.value = lead1_index+1
  lead1.textContent = month_array.at(lead1_index)
  var lead2 = document.querySelector("#lead2")
  lead2.value = lead2_index+1
  lead2.textContent = month_array.at(lead2_index)
  var lead3 = document.querySelector("#lead3")
  lead3.value = lead3_index+1
  lead3.textContent = month_array.at(lead3_index)

  stdate_valor = lead1_index+1
  str_stdate = String(stdate_valor).padStart(2, '0');

  const plot1 = document.querySelector("#plot1") 
  plot1.src = plotsBaseURL+'3-Common_plots/Score-card_'+score_valor+'_'+season_valor+'_3m.png'
  const plot2 = document.querySelector("#plot2") 
  plot2.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.png'
  const plot3 = document.querySelector("#plot3") 
  plot3.src = plotsBaseURL+'2b-EOFs_calc/plots/EOFs/correlations/ERA5_EOF'+eof_num+'_'+season_valor+'.png'  
  const plot4 = document.querySelector("#plot4") 
  plot4.src = plotsBaseURL+folder+folder2+'plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_'+freq+'.'+imgbasename+'.'+season_valor+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
})