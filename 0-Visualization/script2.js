var plotsBaseURL = '../3-Var_modes/'
var imgbasename = '3m.DJF'

var origin = document.querySelector("#sel_origin");
var stdate = document.querySelector("#sel_stdate");
var varn = document.querySelector("#sel_varn");
var score = document.querySelector("#sel_score");

var origin_valor = "ecmwf.s51"
var stdate_valor = "9"
var varn_valor = "nao_box"
var score_valor = "corr"


origin.addEventListener('input',function(evento){
    origin_valor = evento.target.value
    console.log(origin_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
        var folder = "1-Box_calc/"
      } else {
        var folder = "2-EOFs_calc/"
      }
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+folder+'data/plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.png'
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+'data/plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
    })
stdate.addEventListener('input',function(evento){
    stdate_valor = evento.target.value
    console.log(stdate_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
        var folder = "1-Box_calc/"
      } else {
        var folder = "2-EOFs_calc/"
      }
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+folder+'data/plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.png'
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+'data/plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
  })
varn.addEventListener('input',function(evento){
    varn_valor = evento.target.value
    console.log(varn_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
        var folder = "1-Box_calc/"
        var eof_num = "1"
      } else if (varn_valor === "nao_eof") {
        var folder = "2-EOFs_calc/"
        var eof_num = "1"
    } else if (varn_valor === "ea") {
        var folder = "2-EOFs_calc/"
        var eof_num = "2"
    } else if (varn_valor === "eawr") {
        var folder = "2-EOFs_calc/"
        var eof_num = "3"
    } else if (varn_valor === "sca") {
        var folder = "2-EOFs_calc/"
        var eof_num = "4"
    }

    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+folder+'data/plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.png'
    const plot3 = document.querySelector("#plot3") 
    plot3.src = plotsBaseURL+'2-EOFs_calc/data/plots/ERA5_EOF'+eof_num+'.png'  
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+'data/plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
})
score.addEventListener('input',function(evento){
    score_valor = evento.target.value
    console.log(score_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    if (varn_valor === "nao_box") {
      var folder = "1-Box_calc/"
    } else {
      var folder = "2-EOFs_calc/"
    }
  
    const plot1 = document.querySelector("#plot1") 
    plot1.src = plotsBaseURL+'3-Common_plots/Score-card_'+score_valor+'_DJF_3m.png'
    const plot4 = document.querySelector("#plot4") 
    plot4.src = plotsBaseURL+folder+'data/plots/stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'-bootstrap.png'
})
