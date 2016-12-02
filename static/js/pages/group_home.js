$('#group-network h4').click(function() {
    if( $('#network-row').is(":visible")  ) {
      $('#group-network h4 i').removeClass('fa-caret-up').addClass('fa-caret-down');
    } else {
      $('#group-network h4 i').removeClass('fa-caret-down').addClass('fa-caret-up');
    };

    $('#network-row').toggle(function() {
      s.refresh();
      s.refresh();
    });
});


var s = new sigma({
  container: 'network',
  settings: {

    defaultLabelColor: "#000",
    defaultLabelSize: 12,
    defaultLabelBGColor: "#ddd",
    defaultHoverLabelBGColor: "#002147",
    defaultLabelHoverColor: "#fff",

    edgeColor: 'default',
    defaultEdgeColor: '#e7e7e7',

    labelThreshold: 10,
    defaultEdgeType: 'curve',

    hoverFontStyle: "bold",
    fontStyle: "regular",
    activeFontStyle: "regular",

    minNodeSize: 4,
    maxNodeSize: 8,
    minEdgeSize: 1,
    maxEdgeSize: 1,

    zoomMin: .00001,
    zoomMax: 1,
  }
});


var filter = new sigma.plugins.filter(s);

$.ajax({
  'type': 'GET',
  'url': '/api/network/'+ pk +'/',
  'success': function(graph) {
    s.graph.clear();
    s.graph.read(graph);

    s.refresh();
    maxDegree = 0;
    s.graph.nodes().forEach(function(n) {
      maxDegree = Math.max(maxDegree, s.graph.degree(n.id));
    });
    $('#min-degree').attr('max', maxDegree/3)

  },
  'error': function(d) {
    $('#group-network').html("<div class='row'><div class='col-xs-12'><h4 class='text-xs-center'>Network Unavailable <i class='fa fa-exclamation-triangle fa-1'></i></h4></div></div>");
  }
});

s.refresh();
s.bind('hoverNode clickNode', function(e) {
//  console.log(e.type, e.data.node.label, e.data.captor);
  if(e.type == 'clickNode') {
    var node_color = e.data.node.color;
    var second_search = '';
    if (node_color == '#B1FFA8') {
        second_search = 'gene';
    } else if (node_color == '#d1f3ff') {
        second_search = 'disease';
    } else if (node_color == '#ffd1dc') {
        second_search = 'drug';
    } else {
        second_search = '';
    };
    window.open('https://www.google.com/#safe=off&q='+e.data.node.label+'+'+second_search,'_blank');
  }
});

$('#min-degree').on('change', function() {
  var min_degree_num = $(this).val();

  $('#min-degree-val').html(min_degree_num);
  filter
    .undo('min-degree')
    .nodesBy(function(n) {
      return this.degree(n.id) >= min_degree_num;
    }, 'min-degree').apply();
});

$('#reset-btn').on('click', function() {
  filter.undo().apply();
});

$('#network-row i.fa-plus-circle').click(function() {
  var c = s.camera;
  sigma.misc.animation.camera(c, {
    ratio: c.ratio / c.settings('zoomingRatio')
  }, { duration: 250 });
});

$('#network-row i.fa-minus-circle').click(function() {
  var c = s.camera;
  sigma.misc.animation.camera(c, {
    ratio: c.ratio * c.settings('zoomingRatio')
  }, { duration: 250 });
});

$('#network-row i.fa-rotate-right').click(function() {
  var c = s.camera;
  c.goTo({
    angle: c.angle -= .1
  });
});


var draw_dashboard = function(group, quests) {
  $('#group-'+ group.pk).html('');
  var canvas = d3.select('#group-'+ group.pk);

  var available_quests = _.filter(quests, function(item) { return item.enabled && !item.completed });
  var completion_size = _.map(available_quests, function(item) { return item.completions; });

  var completion_scale = d3.scale.linear()
    .domain([_.min(completion_size), _.max(completion_size)])
    .range(['#00CCFF', '#E64C66']);

  var template = _.template( $('#quest-icon-template').html() );
  var attrs = {
    'class': 'quest col-xs-4 col-sm-3 col-md-3 col-lg-2',
  };
  var styles = {
  };
  var quest = canvas.selectAll('.quest').remove();
  var quest = canvas.selectAll('.quest').data(quests);

  quest.enter().append('div')
    .attr(attrs)
    .style(styles)
    .html(function(d, i) {
      return template({
        'd': d,
        'progress': (d.progress.current/d.progress.required)*100,
        'completions_scale': completion_scale(d.completions),
      });
    });
  quest.transition().attr(attrs);
  quest.exit().remove();
};


var group_template = _.template( $('#group-template').html() );
$('#group-selection').append(group_template({'pk': pk}));

$.ajax({
  'type': 'GET',
  'url': '/api/quest/'+ pk +'/',
  'success': function(data) {
    draw_dashboard({'pk': pk}, data);

    $('#group-selection .quest').click(function(evt) {
      var link = $(this).find('a');
      if(link.length) { location.href=link.attr('href'); }
    });

  }
});
