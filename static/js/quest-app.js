var $el = document.querySelector('#score');
od = new Odometer({
  el: $el,
  value: $el.innerHTML,
  format: '(,ddd)',
  theme: 'minimal'
});

var update_score = function() {
    var ajax_settings = {
      url: '/u/points/',
      type: 'GET',
      dataType: 'json',
      success: function(data) {
        od.update(data.points);
      }
    };
    $.ajax(ajax_settings);
};
