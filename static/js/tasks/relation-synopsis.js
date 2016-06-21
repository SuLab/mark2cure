
var document_pk = relation_task_settings.document_pk;
var doc_api_data = null;

$( document ).ready(function(document_pk) {

    $.getJSON('/task/relation/'+ relation_task_settings.document_pk +'/analysis/', function(document_api_data) {
        doc_api_data = document_api_data;

        /* color circles */
        function color_the_circles(doc_api_data) {
            circle_dict = {};
            for(var i = 0; i < doc_api_data.length; i++) {
                var answer_counts;
                var obj = doc_api_data[i];
                var relation_specific_answers = _.filter(doc_api_data, {id: obj.id});
                var answers = relation_specific_answers[0]['answers']
                var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );
                var user_answer = _.filter(answers, {self: true });

                answer_count_nums = _.sortBy(answer_counts, function(d) { return -d['value'] });
                var sum = _.reduce(answer_count_nums, function(memo, num){ return memo + num; }, 0);

                var arr = Object.keys( answer_counts ).map(function ( key ) { return answer_counts[key] });
                var majority_total_votes = Math.max.apply( null, arr );
                majority_answer_id = Object.keys(answer_counts).filter(function(x){ return answer_counts[x] == majority_total_votes; })[0];
                user_answer_id = user_answer[0]['answer']['id']

                percent_of_total = majority_total_votes/sum;

                /* color circles yellow unless the following applies */
                circle = 'yellow';

                /* you match majority and enough data to give green */
                if (user_answer_id == majority_answer_id && percent_of_total >= 0.51){
                    circle = 'green';
                /* you do not match majority and majority did pretty well */
                } else if ( user_answer_id != majority_answer_id && percent_of_total >= 0.51 ) {
                    circle = 'red';
                }
                circle_dict[obj['id']] = circle;

                $('#relation-circle-' + obj['id']).css('color', circle_dict[obj['id']]);
                $('#c1-word-' + obj['id'] ).html("Concept 1: " + obj["concept_a"]["text"]);
                $('#c2-word-' + obj['id'] ).html("Concept 2: " + obj["concept_b"]["text"]);

             }
        };
        circle_dict = color_the_circles(doc_api_data);

        // Click method option
        $( ".relation-circle" ).hover(function () {
            this.style.cursor="pointer";
        });

        var current_id = "";
        var previous_id = "";
        var previous_circle = "";
        var current_circle = "";
        var click_flag = false;

        $( ".relation-circle" ).mousedown(function () {

            previous_id = current_id;
            previous_circle = current_circle;
            if (click_flag == true) {
                previous_circle.style.background="transparent";
            };

            this.style.background="#f7f7f9";
            this.style.cursor="pointer";
            var get_id = (this.id).split("-")[2];
            current_id = get_id;
            relation_pk = current_id;
            function show_results(document_pk, relation_pk, doc_api_data) {
                var relation_specific_answers = _.filter(doc_api_data, {id: parseInt(relation_pk)})
                var answers = relation_specific_answers[0]['answers']
                var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );
                var answer_text = {};
                var personal_ann = '';

                _.each(answers, function(a) {
                  answer_text[a.answer.id] = a.answer.text;
                  if(a['self']) { personal_ann = a.answer.id; }
                });

                var data = [];
                _.each(_.keys(answer_counts), function(answer_key) {
                  data.push({
                    'id': answer_key,
                    'value': answer_counts[answer_key],
                    'label': answer_text[answer_key],
                    'self': answer_key == personal_ann
                  });
                });

                data = _.sortBy(data, function(d) { return -d['value'] });

                var max = d3.sum(data.map(function(i){ return i['value']; }));
                var color = d3.scale.category20c();

                var chart = d3.select('#chart-' + relation_pk).style('width', '100%');
                var bar = chart.selectAll('div')
                  .data(data)
                  .enter()
                    .append('div')
                      .attr('class', 'bar-component')
                      .style('width', function(d) { return ((d['value']/max)*100) + '%'; } )
                      .style('background-color', function(d, i) { return color(i); });

                var bold = function(d, i) {
                  if( i == 1 && d['self'] ) {
                    return '<strong>';
                  } else if( i == 2 && d['self'] ) {
                    return '</strong>';
                  } else {
                    return '';
                  }
                }

                var list = d3.select('#chart-list-' + relation_pk);
                var bar = list.selectAll('li')
                  .data(data)
                  .enter()
                    .append('li')
                      .html(function(d, i) { return '<div class="box" style="background-color:'+color(i)+';"></div> <p>' + bold(d, 1) + ((d['value']/max)*100).toFixed() + '% – ' + d['label'] + bold(d, 2) + '</p>'; });
            }

            show_results(document_pk, relation_pk, doc_api_data);

            current_circle = this;
            $('#relation-' + previous_id).hide();
            $('#relation-' + current_id).show();
            click_flag = true;
        })
    });
});





