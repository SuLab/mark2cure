var NETWORK = {}

NETWORK['init'] = function($el) {
  var width = 600,
      height = 600;

  NETWORK['color'] = d3.scale.category20();
  NETWORK['force'] = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);

  NETWORK['svg'] = d3.select( $el[0] ).append("svg")
    .attr("width", width)
    .attr("height", height);

    d3.selectAll('div.infotip').remove();

    NETWORK['infotip'] = d3.select('body')
      .append('div')
        .attr('class', 'infotip')
        .style('position', 'absolute')
        .style('z-index', '10')
        .style('visibility', 'hidden');

    d3.json("/analysis/network", function(error, graph) { NETWORK.drawNetwork(graph) });

}

NETWORK['drawNetwork'] = function(graph) {
  var self = this;

  NETWORK['force'].stop()
  NETWORK['force']
    .nodes(graph.nodes)
    .links(graph.links);

  var link = NETWORK['svg'].selectAll(".link")
        .data(graph.links)
      .enter().append("line")
        .attr("class", "link")
        .style("stroke-width", function(d) { return Math.sqrt(d.value); });

  var node = NETWORK['svg'].selectAll(".node")
      .data(graph.nodes)
    .enter().append("circle")
      .attr("class", "node")
      .attr("r", 5)
      .style("fill", function(d) { return NETWORK.color(d.group); })
      .call(NETWORK['force'].drag)
      //-- Mouse interactions
      .on('mouseover', function() { return NETWORK['infotip'].style('visibility', 'visible'); })
      .on('mousemove', function(obj) {
        return NETWORK['infotip']
                .style('top', (event.pageY-16)+'px')
                .style('left', (event.pageX+10)+'px')
                .html(function() { return obj.name; });
      })
      .on('mouseout', function() { return NETWORK['infotip'].style('visibility', 'hidden'); });

    node.append("title")
        .text(function(d) { return d.name; });

    NETWORK['force'].on("tick", function() {
      link.attr("x1", function(d) { return d.source.x; })
          .attr("y1", function(d) { return d.source.y; })
          .attr("x2", function(d) { return d.target.x; })
          .attr("y2", function(d) { return d.target.y; });

      node.attr("cx", function(d) { return d.x; })
          .attr("cy", function(d) { return d.y; });
    });

  NETWORK['force'].start()
}

$(document).ready(function() {

  $('tr.annotation-finder').mouseenter(function() {
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


  $('.network').click(function(evt) {
    evt.preventDefault();
    $el = $('#analytics-network')
    $el.show();

    NETWORK.init( $el );
  });


  $('.submit-on-focusout').blur(function(evt) {
    console.log('ran')
    // $(this).closest('form').submit();
  });

  $('.selectize-relationshiptypes').selectize({
      valueField: 'id',
      labelField: 'full_name',
      searchField: ['full_name', 'type', 'id'],
      options: [],
      create: false,
      render: {
          option: function(item, escape) {
              return '<div class="row">' +
                        '<div class="col-md-11 col-md-offset-1">' +
                          '<div class="title">' + escape(item.full_name) + '</div>' +
                        '</div>' +
                     '</div>';
          }
      },
      load: function(query, callback) {
          if (!query.length) return callback();
          $.ajax({
              url: '/document/relationshiptypes/?format=json',
              type: 'GET',
              data: {
                  prefix: query,
                  page_limit: 10
              },
              error: function() {
                  callback();
              },
              success: function(res) {
                  callback(res);
              }
          });
      }
  });




});


// if( this.options.explain ) {
//   self.ui.close.popover({ title   : '',
//                           content : templates.snippets.network_info({}),
//                           html    : true,
//                           trigger : 'manual',
//                           placement : 'down',
//                           container : 'body' });
//   self.ui.close.popover('show');
//   $('button.close-network-popover').click(function(e) {
//     self.ui.close.popover('hide');
//   });

//   //-- Reposition the popover
//   $('.popover').css({'position': 'absolute', top: 120, left: ($('body').width()/2)+100 });
//   this.options.explain = false;
// }
