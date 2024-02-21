var plotsBaseURL = './Score_files/'
var imgbasename = 'DJF_3m'

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
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_'+varn_valor+'_'+imgbasename+'.png'
})
stdate.addEventListener('input',function(evento){
    stdate_valor = evento.target.value
    console.log(stdate_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_'+varn_valor+'_'+imgbasename+'.png'
})
varn.addEventListener('input',function(evento){
    varn_valor = evento.target.value
    console.log(varn_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');
    
    const plot2 = document.querySelector("#plot2") 
    plot2.src = plotsBaseURL+origin_valor.split('.')[0]+'_'+origin_valor.split('.')[1]+'_stmonth'+str_stdate+'_'+varn_valor+'_'+imgbasename+'.png'
    const plot3 = document.querySelector("#plot3") 
    plot3.src = plotsBaseURL+'ERA5_'+varn_valor+'.png'  
})
score.addEventListener('input',function(evento){
    score_valor = evento.target.value
    console.log(score_valor)
    var str_stdate = String(stdate_valor).padStart(2, '0');
    
    const plot1 = document.querySelector("#plot1") 
    plot1.src = plotsBaseURL+score_valor+'_'+imgbasename+'.png'
})
