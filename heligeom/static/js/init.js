(function($){
  $(function(){

    $('.sidenav').sidenav();

  }); // end of document ready
})(jQuery); // end of jQuery name space

$(document).ready(function(){
  $('.collapsible').collapsible();
});

// When called, remove the "hide" class of one element id given in parameter
function load(ElementId){
  document.getElementById(ElementId).classList.remove("hide");
};

function showErrorLogs() {
  var text = document.getElementById("logs");
  if (text.style.display === "none") {
    text.style.display = "block";
  } else {
    text.style.display = "none";
  }
}


// Limit file input size
function checkSize(){
  if(this.files[0].size > 10097152){
    alert("Maximum size accepted is 10Mb.");
    this.value = "";
 };
}

var uploadField = document.getElementById("input_file");
uploadField.onchange = checkSize;
var uploadField2 = document.getElementById("input_file_2nd");
uploadField2.onchange = checkSize;



// Pre populate form with JS
document.getElementById("example").onclick = function() {
  //Fill PDBid
  document.getElementById("pdb_id").value="2GLS";
  //Fill Chains
  document.getElementById("chain1_id").value="A";
  document.getElementById("chain2_id").value="B";
}

// Use to create a "loading page"
function loading(){
    $("#loading").show();
    $("#content_run").hide();
}

