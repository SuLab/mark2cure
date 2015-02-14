var drawUserWithGolden = function() {
  var $sections = $('p.paragraph');

  $.each($sections, function() {
    var section = $(this);
    var $words = section.find('span');

    $.each($words, function(w) {
      var word = $(this),
          next_word = word.next(),

          gm_ann = word.data('gmannid') || 'None',
          next_word_gm_ann = next_word.data('gmannid') || 'None',

          u_ann = word.data('uannid') || 'None',
          next_word_u_ann = next_word.data('uannid') || 'None';

      if(gm_ann != 'None') { word.addClass('golden'); }
      if(u_ann != 'None') { word.addClass('user_annotated'); }

      if( gm_ann != 'None' && gm_ann != next_word_gm_ann ) {
        word.html( word.html().trim() );
        word.addClass('neighbor_gm');
        word.css({'margin-right': '5px'});

        /* Check if it messes up user anns */
        if( u_ann != 'None' && u_ann == next_word_u_ann) {
          var $putty = $("<span class='putty-rect user_annotated'></span>").appendTo(section),
              section_offset = section.offset(),
              word_offset = word.offset(),
              word_pos = word.position();
          $putty.css({'top' : word_pos.top, 'left' : word_pos.left + word.width() - 1, 'height': word.height() });
        }
      };

      if( u_ann != 'None' && u_ann != next_word_u_ann ) {
        word.html( word.html().trim() );
        word.addClass('neighbor_u');
        word.css({'margin-right': '5px'});

        /* Check if it messes up GM anns */
        if( gm_ann != 'None' && gm_ann == next_word_gm_ann) {
          var $putty = $("<span class='putty-rect golden'></span>").appendTo(section),
              position = word.position();
          $putty.css({'top' : position.top + word.height(), 'left' : position.left + word.width() - 1, 'height': '6px' });
        }

      }
    });

  });
};

$(document).ready(function() {
  $('#simple-instructions-open').click(function() {
    $('#accordion').slideDown();
  });

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

  $('.enter-form-submit input[type=password]').keypress(function(evt) {
    if(evt.which == 13) {
      evt.preventDefault();
      evt.stopPropagation();
      $form = $(this).closest('.modal-content').find('form');
      $form.submit();
    }
  });

  $('.modal-form-submit').on('click', function(evt) {
    evt.preventDefault();
    evt.stopPropagation();
    $form = $(this).closest('.modal-content').find('form');
    $form.submit();
  });

  $('.form-via-ajax').on('submit', function(evt) {
    evt.preventDefault();
    evt.stopPropagation();

    var form = $(this),
        data = {},
        input;

    var $textfield = form.find("textarea[name='text']");
    if($textfield.length) {
      if( $textfield.val().length == 0 ){
        $textfield.parent().addClass('has-error');
        return;
      };
    }

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
