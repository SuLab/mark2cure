define(['marionette', 'templates', 'vent',
        'd3'],

        function (Marionette, templates, vent) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.analytics.network,
    templateHelpers : function() { return this.options; },
    className : 'analytics-view',

    ui : {
      network : '#network',
      graph   : 'svg',
      close   : 'button.close',
      refresh : 'i.refresh'
    },

    events : {
      'click button.close'  : 'closeAnalytics',
      'click i.refresh'     : 'refreshGraph'
    },

    initialize : function() {
      var self = this;
      vent.on('network:refresh', function() {
        self.refreshGraph();
      });
    },

    onRender : function() {
      var self = this,
          width = 600,
          height = 260;

      if( this.options.explain ) {
        self.ui.close.popover({ title   : '',
                                content : templates.snippets.network_info({}),
                                html    : true,
                                trigger : 'manual',
                                placement : 'down',
                                container : 'body' });
        self.ui.close.popover('show');
        $('button.close-network-popover').click(function(e) {
          self.ui.close.popover('hide');
        });

        //-- Reposition the popover
        $('.popover').css({'position': 'absolute', top: 120, left: ($('body').width()/2)+100 });
        this.options.explain = false;
      }

      this.options.network = {};
      this.options.network.color = d3.scale.category20();
      this.options.network.force = d3.layout.force()
          .charge(-120)
          .linkDistance(30)
          .size([width, height]);

      this.options.network.svg = d3.select( this.ui.network[0] ).append("svg")
          .attr("width", width)
          .attr("height", height);

      d3.selectAll('div.infotip').remove();

      this.options.network.infotip = d3.select('body')
        .append('div')
          .attr('class', 'infotip')
          .style('position', 'absolute')
          .style('z-index', '10')
          .style('visibility', 'hidden');

      d3.json("/api/v1/network", function(error, graph) { self.drawNetwork(graph) });
    },

    //
    //-- Events
    //
    refreshGraph : function(evt) {
      var self = this;
      // d3.json("/api/v1/network", function(error, graph) { self.drawNetwork(graph) });
      this.render();
    },

    closeAnalytics : function(evt) {
      evt.preventDefault();
      vent.trigger('navigate:analytics', {toggle: false});
    },

    drawNetwork : function(graph) {
      var self = this;
      this.options.network.force.stop()
      this.options.network.force
        .nodes(graph.nodes)
        .links(graph.links);

      var link = this.options.network.svg.selectAll(".link")
            .data(graph.links)
          .enter().append("line")
            .attr("class", "link")
            .style("stroke-width", function(d) { return Math.sqrt(d.value); });

      var node = this.options.network.svg.selectAll(".node")
          .data(graph.nodes)
        .enter().append("circle")
          .attr("class", "node")
          .attr("r", 5)
          .style("fill", function(d) { return self.options.network.color(d.group); })
          .call(self.options.network.force.drag)
          //-- Mouse interactions
          .on('mouseover', function() { return self.options.network.infotip.style('visibility', 'visible'); })
          .on('mousemove', function(obj) {
            return self.options.network.infotip
                    .style('top', (event.pageY-16)+'px')
                    .style('left', (event.pageX+10)+'px')
                    .html(function() { return obj.name; });
          })
          .on('mouseout', function() { return self.options.network.infotip.style('visibility', 'hidden'); });

        node.append("title")
            .text(function(d) { return d.name; });

        this.options.network.force.on("tick", function() {
          link.attr("x1", function(d) { return d.source.x; })
              .attr("y1", function(d) { return d.source.y; })
              .attr("x2", function(d) { return d.target.x; })
              .attr("y2", function(d) { return d.target.y; });

          node.attr("cx", function(d) { return d.x; })
              .attr("cy", function(d) { return d.y; });
        });

      this.options.network.force.start()
    }

  });
});
