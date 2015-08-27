$(document).ready(function() {
  var $el = document.querySelector('#score');
  if($el) {
    od = new Odometer({
      el: $el,
      value: $el.innerHTML,
      format: '(,ddd)',
      theme: 'minimal'
    });
  }
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

var drawUserFScoreLine = function(div_sel_str, group_pk) {
    var margin = {top: 30, right: 40, bottom: 30, left: 50},
        width = 500 - margin.left - margin.right,
        height = 270 - margin.top - margin.bottom;

    var parseDate = d3.time.format("%Y-%m-%dT%H:%M:%SZ").parse;

    var x = d3.time.scale().range([0, width]);
    var y0 = d3.scale.linear().range([height, 0]);
    var y1 = d3.scale.linear().range([height, 0]);

    var xAxis = d3.svg.axis().scale(x)
      .orient("bottom").ticks(5);

    var yAxisLeft = d3.svg.axis().scale(y0)
      .orient("left").ticks(5);

    var yAxisRight = d3.svg.axis().scale(y1)
      .orient("right").ticks(5);

    var valueline = d3.svg.line()
      .x(function(d) { return x(d.date); })
      .y(function(d) { return y0(d.fscore); });

    var valueline2 = d3.svg.line()
      .x(function(d) { return x(d.date); })
      .y(function(d) { return y1(d.pairings); });

    var svg = d3.select(div_sel_str)
      .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    d3.json('/api/analysis/group/'+group_pk+'/user/?format=json', function(error, data) {

      data.forEach(function(d) {
        d.date = parseDate(d.created);
        d.fscore = d['f-score'];
      });


      x.domain(d3.extent(data, function(d) { return d.date; }));
      y0.domain([0, 1]);
      y1.domain([0, _.max(_.pluck(data, 'pairings'))*1.2 ]);

      svg.append("path")
        .style("stroke", "#7F3CFF")
        .attr("d", valueline(data));

      svg.append("path")
        .style("stroke", "#E85997")
        .attr("d", valueline2(data));

      svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

      svg.append("g")
        .attr("class", "y axis")
        .style("fill", "#7F3CFF")
        .call(yAxisLeft);

      svg.append("g")
        .attr("class", "y axis")
        .attr("transform", "translate(" + width + " ,0)")
        .style("fill", "#E85997")
        .call(yAxisRight);
    });

}

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
