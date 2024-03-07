var plotsBaseURL = '../1-Sf_variables/data/plots/'
var imgbasename = '3m.DJF'

var origin = document.querySelector("#sel_origin");
var stdate = document.querySelector("#sel_stdate");
var varn = document.querySelector("#sel_varn");
var score = document.querySelector("#sel_score");

var origin_valor = "ecmwf.s51"
var stdate_valor = "11"
var varn_valor = "t2m"
var score_valor = "corr"


origin.addEventListener('input',function(evento){
    origin_valor = evento.target.value
    console.log(origin_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+'stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'.png'
})
stdate.addEventListener('input',function(evento){
    stdate_valor = evento.target.value
    console.log(stdate_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+'stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'.png'
})
varn.addEventListener('input',function(evento){
    varn_valor = evento.target.value
    console.log(varn_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+'stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'.png'
})
score.addEventListener('input',function(evento){
    score_valor = evento.target.value
    console.log(score_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');

    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+'stmonth'+str_stdate+'/'+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_hindcast1993-2016_monthly.'+imgbasename+'.'+varn_valor+'.'+score_valor+'.png'
    const plot1 = document.querySelector("#plot1") 
    plot1.src = plotsBaseURL+'Score-card_'+score_valor+'_DJF_3m.png'
})
