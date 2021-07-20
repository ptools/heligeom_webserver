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
