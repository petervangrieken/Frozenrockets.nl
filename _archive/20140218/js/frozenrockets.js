
window.addEventListener('devicelight', function(e) {
  var lux = e.value;

  if(lux < 50) { // dim
    document.body.className = 'dim';
  }
  if(lux >= 50 && lux <= 1000) {
    document.body.className = 'normal';
  }
  if(lux > 1000)  { // bright
    document.body.className = 'bright';
  } 
});



///if (window.DeviceMotionEvent!=undefined) {

        var angle = 0;
        var multiplier = 9 / 90;
        window.ondeviceorientation = function(event) {
            var new_angle = ((self.orientation === "portrait" ) ? event.beta : event.gamma) * multiplier;
            if(new_angle != angle){
              $('#angle').text(new_angle);
                angle = new_angle;
                $('.home-content .home-background').css({
                    transform: 'translate3d(0, '+angle+'% , 0)'
                })
                $('.home-content .home-link-container .js-home-link .buttonlink').css({
                    transform: 'translate3d(0, '+(angle*-1)+'px , 0)'
                });
            }
        }

//}