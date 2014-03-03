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