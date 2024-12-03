M.AutoInit();

(function($){
  $(function(){

    $('.sidenav').sidenav();

  }); // end of document ready
})(jQuery); // end of jQuery name space

$(document).ready(function(){
  $('.collapsible').collapsible();
  $('.tooltipped').tooltip();
  $('.scrollspy').scrollSpy();
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




// Use to create a "loading page"
function loading(){
    $("#loading").show();
    $("#content_run").hide();
}

