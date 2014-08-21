$(document).ready(function() {

  $('.annotation-finder').mouseenter(function() {
    var $ann = $(this),
        $section = $('#'+ $ann.data('section')+'.paragraph')
        needle_start = $ann.data('start'),
        needle_length = $ann.data('text').length,
        start = 0;

    $.each( $section.find('span'), function() {
      start = $(this).data('starti');

      if(start >= needle_start && start <= needle_start+needle_length) {
        $(this).addClass('focused');
      }

    });
  }).mouseleave(function() {
    var $section = $('#'+ $(this).data('section')+'.paragraph');
    $.each( $section.find('span.focused'), function() {
      $(this).removeClass('focused');
    });
  });

  $('.form-via-ajax').on('submit', function(evt) {
    evt.preventDefault();
    evt.stopPropagation();

    var form = $(this),
        data = {},
        input;

    $.each( form.find('input, textarea, select'), function() {
      input = $(this);
      data[ input.attr('name') ] = input.val();
      input.val('')
    });

    $.ajax({
      type: 'POST',
      url: form.attr('action'),
      data: data,
      cache: false,
      async: false,
      success: function() {
        if( form.parent().hasClass('modal-body') ) {
          form.closest('.modal').modal('hide');
        }
      }
    });
  });

  $('.submit-on-focusout').blur(function(evt) {
    console.log('ran')
    // $(this).closest('form').submit();
  });

  $('.selectize-basic').selectize();
  $('.selectize-occupation').selectize({
    valueField: 'name',
    labelField: 'name',
    searchField: 'name',
    options: [
      {name: 'Student'},
      {name: 'Science'},
      {name: 'Biological Sciences'},
      {name: 'Chemical Sciences'},
      {name: 'Computer'},
      {name: 'Labor'},
      {name: 'Unemployed'},
      {name: 'Technical'},
      {name: 'Programmer'},
      {name: 'Business'},
      {name: 'Finance'},
      {name: 'Legal'},
      {name: 'Art'},
      {name: 'Education'},
      {name: 'Retired'},
    ],
    create: true,
    maxItems: 5,
  })

  $('.selectize-motivation').selectize({
    valueField: 'name',
    labelField: 'name',
    searchField: 'name',
    options: [
      {name: 'I want to make money'},
      {name: 'I want entertainment'},
      {name: 'I want to help science'},
    ],
    create: true,
    maxItems: 5,
  })

});
