// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, true, false ];
var arrayMetadata    = [ [ "1", "GSM1625995", "con6h-rep1", "infection: control", "time: 6 hpi", "Gene expression after 6h in control-infected U251", "control", "6 hpi" ], [ "2", "GSM1625996", "con6h-rep2", "infection: control", "time: 6 hpi", "Gene expression after 6h in control-infected U251", "control", "6 hpi" ], [ "3", "GSM1625997", "con6h-rep3", "infection: control", "time: 6 hpi", "Gene expression after 6h in control-infected U251", "control", "6 hpi" ], [ "4", "GSM1625998", "con12h-rep1", "infection: control", "time: 12 hpi", "Gene expression after 12h in control-infected U251", "control", "12 hpi" ], [ "5", "GSM1625999", "con12h-rep2", "infection: control", "time: 12 hpi", "Gene expression after 12h in control-infected U251", "control", "12 hpi" ], [ "6", "GSM1626000", "con12h-rep3", "infection: control", "time: 12 hpi", "Gene expression after 12h in control-infected U251", "control", "12 hpi" ], [ "7", "GSM1626001", "con24h-rep1", "infection: control", "time: 24 hpi", "Gene expression after 24h in control-infected U251", "control", "24 hpi" ], [ "8", "GSM1626002", "con24h-rep2", "infection: control", "time: 24 hpi", "Gene expression after 24h in control-infected U251", "control", "24 hpi" ], [ "9", "GSM1626003", "con24h-rep3", "infection: control", "time: 24 hpi", "Gene expression after 24h in control-infected U251", "control", "24 hpi" ], [ "10", "GSM1626004", "hm6h-rep1", "infection: H5N1", "time: 6 hpi", "Gene expression after 6h in H5N1-infected U251", "H5N1", "6 hpi" ], [ "11", "GSM1626005", "hm6h-rep2", "infection: H5N1", "time: 6 hpi", "Gene expression after 6h in H5N1-infected U251", "H5N1", "6 hpi" ], [ "12", "GSM1626006", "hm6h-rep3", "infection: H5N1", "time: 6 hpi", "Gene expression after 6h in H5N1-infected U251", "H5N1", "6 hpi" ], [ "13", "GSM1626007", "hm12h-rep1", "infection: H5N1", "time: 12 hpi", "Gene expression after 12h in H5N1-infected U251", "H5N1", "12 hpi" ], [ "14", "GSM1626008", "hm12h-rep2", "infection: H5N1", "time: 12 hpi", "Gene expression after 12h in H5N1-infected U251", "H5N1", "12 hpi" ], [ "15", "GSM1626009", "hm12h-rep3", "infection: H5N1", "time: 12 hpi", "Gene expression after 12h in H5N1-infected U251", "H5N1", "12 hpi" ], [ "16", "GSM1626010", "hm24h-rep1", "infection: H5N1", "time: 24 hpi", "Gene expression after 24h in H5N1-infected U251", "H5N1", "24 hpi" ], [ "17", "GSM1626011", "hm24h-rep2", "infection: H5N1", "time: 24 hpi", "Gene expression after 24h in H5N1-infected U251", "H5N1", "24 hpi" ], [ "18", "GSM1626012", "hm24h-rep3", "infection: H5N1", "time: 24 hpi", "Gene expression after 24h in H5N1-infected U251", "H5N1", "24 hpi" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
